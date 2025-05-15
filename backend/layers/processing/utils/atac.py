import os
import scanpy as sc
import pandas as pd
import pysam
import tiledb
import numpy as np
import json
from collections import defaultdict
import logging


logger = logging.getLogger(__name__)


# Config
# FRAGMENT_FILE = "8b7ffdda-ad10-4985-a37e-436a79cce4dd-fragment.tsv.bgz"
# TILEDB_GROUP = "atac_coverage_tiledb_with_cell_id"
# COVERAGE_ARRAY = f"{TILEDB_GROUP}/coverage"
USE_NORMALIZATION = False
MAX_BINS = 100_000
BIN_SIZE = 10000
CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569
}


class ATACDataProcessor:
    """
    ATACDataProcessor is a class designed to process ATAC-seq fragment files and generate binned coverage data for downstream analysis. 
    It provides methods to extract metadata, map cell names to IDs, create a TileDB array schema, and write binned coverage data.
    Attributes:
        fragment_artifact_id (str): The path to the fragment file to be processed.
        ctx (tiledb.Ctx): TileDB context for managing TileDB operations.
        bin_size (int): The size of bins used for coverage calculation.
        chrom_lengths (dict): A dictionary mapping chromosome names to their lengths.
        bins_per_chrom (dict): A dictionary mapping chromosome names to the number of bins.
        scale_factors (dict): A dictionary mapping chromosome names to scale factors for normalization.
        max_bins (int): The maximum number of bins per chromosome.
    Methods:
        __init__(fragment_artifact_id=None, ctx=None):
            Initializes the ATACDataProcessor with optional fragment file path and TileDB context.
        extract_cell_metadata_from_h5ad(obs: pd.DataFrame, obs_column: str = "cell_type") -> pd.DataFrame:
            Extracts cell metadata from the H5AD file and returns a DataFrame with cell names and metadata.
        map_cell_name_to_ids() -> Tuple[int, dict, int, dict]:
            Maps cell names to unique IDs and returns the maximum chromosome ID, cell ID map, number of cells, and chromosome map.
        create_dataframe_array(array_name: str, max_chrom: int, num_cells: int):
            Creates a TileDB sparse array schema for storing binned coverage data.
        write_binned_coverage_per_chrom(array_name: str, chrom_map: dict, cell_id_map: dict):
            Writes binned coverage data for each chromosome to the TileDB array.
        process_fragment_file(obs: pd.DataFrame, array_name: str) -> Tuple[pd.DataFrame, dict]:
            Processes the fragment file, generates binned coverage data, and returns cell metadata and cell ID map.
    """

    def __init__(self, fragment_artifact_id=None, ctx=None):
        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.chrom_lengths = CHROM_LENGTHS
        # ######
        # bins = (248956422 + 10000 - 1) // 10000
        #  = 248966421 // 10000
        #  = 24896
        # ######
        self.bins_per_chrom = {
            chrom: int(length / self.bin_size) + 1 for chrom, length in self.chrom_lengths.items()
        }
        self.scale_factors = {chrom: max(1, (bins + MAX_BINS - 1)// MAX_BINS) if USE_NORMALIZATION else 1 for chrom, bins in self.bins_per_chrom.items()}
        self.max_bins = MAX_BINS if USE_NORMALIZATION else max(self.bins_per_chrom.values()) - 1

    def extract_cell_metadata_from_h5ad(self, obs: pd.DataFrame, obs_column:str = "cell_type"):
        """
        Extract cell metadata from the H5AD file.
        """
        if obs_column not in obs.columns:
            raise ValueError(f"Column {obs_column} not found in obs DataFrame.")
        df = obs[[obs_column]].copy()
        df.index.name = "cell_name"
        return df.reset_index()
    
    def map_cell_name_to_ids(self):
        """
        Get the maximum chromosome ID.
        """
        chrom_map = defaultdict(lambda: 0)
        for i, chrom in enumerate(self.chrom_lengths, 1):
            chrom_map[chrom] = i
        max_chrom = max(chrom_map.values())
    
        tabix = pysam.TabixFile(self.fragment_artifact_id)
        cell_id_map = {}
        next_id = 0
        for chrom_str in chrom_map.keys():
            try:
                for row in tabix.fetch(chrom_str):
                    cell = row.strip().split('\t')[3]
                    if cell not in cell_id_map:
                        cell_id_map[cell] = next_id
                        next_id += 1
            except ValueError:
                continue
        num_cells = len(cell_id_map)

        return max_chrom, cell_id_map, num_cells, chrom_map

    def create_dataframe_array(self, array_name, max_chrom, num_cells):
        domain = tiledb.Domain(
            tiledb.Dim(name="chrom", domain=(1, max_chrom), tile=1, dtype=np.uint32),
            tiledb.Dim(name="bin", domain=(0, self.max_bins), tile=10, dtype=np.uint32),
            tiledb.Dim(name="cell_id", domain=(0, num_cells - 1), tile=1000, dtype=np.uint32),
        )

        schema = tiledb.ArraySchema(
            domain=domain,
            attrs=[tiledb.Attr(name="coverage", dtype=np.float32)],
            sparse=True,
            allows_duplicates=False,
        )
        tiledb.SparseArray.create(array_name, schema)

    def write_binned_coverage_per_chrom(self, array_name, chrom_map, cell_id_map):
        tabix = pysam.TabixFile(self.fragment_artifact_id)
        # chrome_to_ingest = chrom_map.keys()
        chrome_to_ingest = ["chr1"]


        with tiledb.SparseArray(array_name, mode="w", ctx=self.ctx) as A:
            for chrom_str in chrome_to_ingest:
                chrom_id = chrom_map[chrom_str]
                print(f"Processing {chrom_str}...")
                try:
                    rows = list(tabix.fetch(chrom_str))
                except ValueError:
                    continue

                scale_factor = self.scale_factors.get(chrom_str, 1)
                data = []
                for row in rows:
                    _, start, _, cell, _ = row.strip().split('\t')
                    raw_bin = int(start) // self.bin_size
                    norm_bin = raw_bin // scale_factor
                    cell_idx = cell_id_map[cell]
                    data.append((chrom_id, norm_bin, cell_idx))

                if not data:
                    continue
            
                df = pd.DataFrame(data, columns=["chrom", "bin", "cell_id"])
                binned = df.groupby(["chrom", "bin", "cell_id"]).size().reset_index(name="coverage")
                A[(binned["chrom"].values, binned["bin"].values, binned["cell_id"].values)] = {
                    binned["coverage"].values
                }


    def process_fragment_file(self, obs: pd.DataFrame, array_name: str):
        """
        Process the fragment file and return a DataFrame with the processed data.
        """
        df_meta = self.extract_cell_metadata_from_h5ad(obs)
        max_chrom, cell_id_map, num_cells, chrom_map = self.map_cell_name_to_ids()
        self.create_dataframe_array(array_name, max_chrom, num_cells)
        self.write_binned_coverage_per_chrom(array_name, chrom_map, cell_id_map)

        return df_meta, cell_id_map
        