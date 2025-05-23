import os
import scanpy as sc
import pandas as pd
import pysam
import tiledb
import numpy as np
import json
from collections import defaultdict
import logging
import time


logger = logging.getLogger(__name__)


# Config
BIN_SIZE = 100
CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569
}


class ATACDataProcessor:
    def __init__(self, fragment_artifact_id=None, ctx=None):
        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.chrom_lengths = CHROM_LENGTHS

        self.bins_per_chrom = {
            chrom: int(length / self.bin_size) + 1 for chrom, length in self.chrom_lengths.items()
        }
        self.max_bins = max(self.bins_per_chrom.values()) - 1

    def extract_cell_metadata_from_h5ad(self, obs: pd.DataFrame, obs_column: str = "cell_type"):
        if obs_column not in obs.columns:
            raise ValueError(f"Column {obs_column} not found in obs DataFrame.")
        df = obs[[obs_column]].copy()
        df = df.rename_axis("cell_name").reset_index()
        # print("df", df)
        # df.to_csv("output.csv", index=True)
        return df

    def map_cell_type_info(self, valid_barcodes: set, cell_type_map: dict):
        chrom_map = defaultdict(lambda: 0)
        for i, chrom in enumerate(self.chrom_lengths, 1):
            chrom_map[chrom] = i
        max_chrom = max(chrom_map.values())

        tabix = pysam.TabixFile(self.fragment_artifact_id)
        cell_id_map = {}

        for chrom_str in chrom_map.keys():
            try:
                for row in tabix.fetch(chrom_str):
                    # print("tabix row", row)
                    cell = row.strip().split('\t')[3]
                    if cell not in valid_barcodes:
                        continue
                    if cell not in cell_id_map:
                        cell_id_map[cell] = {
                            "cell_type": cell_type_map[cell],
                        }
                        # print("cell_id_map", cell_id_map)
            except ValueError:
                continue

        return max_chrom, cell_id_map, chrom_map


    def create_dataframe_array(self, array_name, max_chrom):
        domain = tiledb.Domain(
            tiledb.Dim(name="chrom", domain=(1, max_chrom), tile=1, dtype=np.uint32),
            tiledb.Dim(name="bin", domain=(0, self.max_bins), tile=10, dtype=np.uint32),
            tiledb.Dim(name="cell_type", dtype="ascii"),  # cell_type as dimension
        )

        schema = tiledb.ArraySchema(
            domain=domain,
        attrs=[
            tiledb.Attr(name="coverage", dtype=np.float32),
            tiledb.Attr(name="total_coverage", dtype=np.float32),
            tiledb.Attr(name="normalized_coverage", dtype=np.float32)
        ],
            sparse=True,
            allows_duplicates=False,
        )
        tiledb.SparseArray.create(array_name, schema)


    def write_binned_coverage_per_chrom(self, array_name, chrom_map, cell_id_map):
        tabix = pysam.TabixFile(self.fragment_artifact_id)
        chrome_to_ingest = list(chrom_map.keys())

        all_binned_data = []

        for chrom_str in chrome_to_ingest:
            chrom_id = chrom_map[chrom_str]
            print(f"Processing {chrom_str}...")
            try:
                rows = list(tabix.fetch(chrom_str))
            except ValueError:
                continue

            data = []
            for row in rows:
                _, start, end, cell, _ = row.strip().split('\t')
                if cell not in cell_id_map:
                    continue
                raw_bin_start = int(start) // self.bin_size
                raw_bin_end = int(end) // self.bin_size

                cell_type = cell_id_map[cell]["cell_type"]
                data.append((chrom_id, raw_bin_start, cell_type))
                data.append((chrom_id, raw_bin_end, cell_type))

            if not data:
                continue

            df = pd.DataFrame(data, columns=["chrom", "bin", "cell_type"])
            binned = df.groupby(["chrom", "bin", "cell_type"]).size().reset_index(name="coverage")
            all_binned_data.append(binned)

        if not all_binned_data:
            return

        # Concatenate all chromosome data
        full_df = pd.concat(all_binned_data, ignore_index=True)
        # full_df.to_csv("full_df.csv", index=True)

        # Compute total coverage across all chromosomes
        cell_type_totals = full_df.groupby("cell_type")["coverage"].sum()
        normalization_factor = 200_000
        full_df["total_coverage"] = full_df["cell_type"].map(cell_type_totals)
        full_df["normalized_coverage"] = (full_df["coverage"] / full_df["total_coverage"]) * normalization_factor

        with tiledb.SparseArray(array_name, mode="w", ctx=self.ctx) as A:
            A[
                (
                    full_df["chrom"].astype("int32").to_numpy(),
                    full_df["bin"].astype("int32").to_numpy(),
                    full_df["cell_type"].to_numpy(),
                )
            ] = {
                "coverage": full_df["coverage"].astype("int32").to_numpy(),
                "total_coverage": full_df["total_coverage"].astype("int32").to_numpy(),
                "normalized_coverage": full_df["normalized_coverage"].to_numpy(),
            }



    def process_fragment_file(self, obs: pd.DataFrame, array_name: str):
        df_meta = self.extract_cell_metadata_from_h5ad(obs)
        valid_barcodes = set(df_meta["cell_name"])
        cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"]))

        max_chrom, cell_id_map, chrom_map = self.map_cell_type_info(valid_barcodes, cell_type_map)

        dropped = len(valid_barcodes - set(cell_id_map.keys()))
        if dropped:
            logger.warning(f"{dropped} barcodes in .h5ad were not found in fragment file.")

        self.create_dataframe_array(array_name, max_chrom)
        self.write_binned_coverage_per_chrom(array_name, chrom_map, cell_id_map)

        return df_meta, cell_id_map


