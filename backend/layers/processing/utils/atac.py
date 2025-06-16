import logging
import os
from collections import defaultdict
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd
import pysam
import tiledb

logger = logging.getLogger(__name__)


# Config
BIN_SIZE = 100
CHROM_LENGTHS = {
    "hg38": {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr21": 46709983,
        "chr22": 50818468,
        "chrX": 156040895,
        "chrY": 57227415,
        "chrM": 16569,
    },
    "mm39": {
        "chr1": 195154279,
        "chr2": 181755017,
        "chr3": 159745316,
        "chr4": 156860686,
        "chr5": 151758149,
        "chr6": 149588044,
        "chr7": 144995196,
        "chr8": 130127694,
        "chr9": 124359700,
        "chr10": 130530862,
        "chr11": 122082543,
        "chr12": 120129022,
        "chr13": 120421639,
        "chr14": 124902244,
        "chr15": 104043685,
        "chr16": 98207768,
        "chr17": 94987271,
        "chr18": 90702639,
        "chr19": 61431566,
        "chrX": 169476879,
        "chrY": 91455967,
        "chrM": 16299,
    },
}
NORMALIZATION_FACTOR = 2_000_000


class ATACDataProcessor:
    def __init__(self, fragment_artifact_id: Optional[str] = None, ctx: Optional[tiledb.Ctx] = None) -> None:
        if fragment_artifact_id is not None and not os.path.exists(fragment_artifact_id):
            raise FileNotFoundError(f"Fragment file not found: {fragment_artifact_id}")

        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.chrom_lengths = CHROM_LENGTHS
        self.normalization_factor = NORMALIZATION_FACTOR

    @staticmethod
    def get_genome_version(organism_ontology_term_id: str) -> str:
        mapping = {
            "NCBITaxon:9606": "hg38",
            "NCBITaxon:10090": "mm39",
        }

        if organism_ontology_term_id not in mapping:
            raise ValueError(f"Unknown organism ontology term ID: {organism_ontology_term_id}")
        return mapping[organism_ontology_term_id]

    def extract_cell_metadata_from_h5ad(
        self, obs: pd.DataFrame, obs_column: str = "cell_type"
    ) -> Tuple[pd.DataFrame, str]:
        if obs_column not in obs.columns:
            raise ValueError(f"Column {obs_column} not found in obs DataFrame.")

        if "organism_ontology_term_id" not in obs.columns:
            raise ValueError("Column 'organism_ontology_term_id' is required but not found in obs DataFrame.")

        organism_ontology_term_id = obs["organism_ontology_term_id"].iloc[0]
        if pd.isna(organism_ontology_term_id):
            raise ValueError("organism_ontology_term_id cannot be null/NaN.")

        df = obs[[obs_column]].copy()
        df = df.rename_axis("cell_name").reset_index()
        genome_version = self.get_genome_version(organism_ontology_term_id)
        return df, genome_version

    def calculate_max_bins(self, genome_version: str) -> int:
        bins_per_chrom = {
            chrom: int(length / self.bin_size) + 1 for chrom, length in self.chrom_lengths[genome_version].items()
        }
        return max(bins_per_chrom.values()) - 1

    def build_chrom_mapping(self, genome_version: str) -> Tuple[int, Dict[str, int]]:
        chrom_map = defaultdict(lambda: 0)
        for i, chrom in enumerate(self.chrom_lengths[genome_version], 1):
            chrom_map[chrom] = i
        max_chrom = max(chrom_map.values())
        return max_chrom, chrom_map

    def create_dataframe_array(self, array_name: str, max_chrom: int, max_bins: int) -> None:
        compression = tiledb.FilterList(
            [
                tiledb.BitShuffleFilter(),
                tiledb.ZstdFilter(level=3),
            ]
        )

        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=3)])

        domain = tiledb.Domain(
            tiledb.Dim(name="chrom", domain=(1, max_chrom), tile=1, dtype=np.uint32, filters=dim_filters),
            tiledb.Dim(name="bin", domain=(0, max_bins), tile=10, dtype=np.uint32, filters=dim_filters),
            tiledb.Dim(name="cell_type", dtype="ascii", filters=dim_filters),
        )

        schema = tiledb.ArraySchema(
            domain=domain,
            attrs=[
                tiledb.Attr(name="coverage", dtype=np.int32, filters=compression),
                tiledb.Attr(name="total_coverage", dtype=np.int32, filters=compression),
                tiledb.Attr(name="normalized_coverage", dtype=np.float32, filters=compression),
            ],
            sparse=True,
            allows_duplicates=False,
        )
        tiledb.SparseArray.create(array_name, schema)

    def write_binned_coverage_per_chrom(
        self, array_name: str, chrom_map: Dict[str, int], cell_type_map: Dict[str, str], valid_barcodes: Set[str]
    ) -> None:
        chrome_to_ingest = list(chrom_map.keys())
        found_cells = set()  # Track cells we've actually found

        # Use incremental aggregation to avoid memory explosion
        coverage_aggregator = defaultdict(int)  # (chrom, bin, cell_type) -> coverage_count

        with pysam.TabixFile(self.fragment_artifact_id) as tabix:
            for chrom_str in chrome_to_ingest:
                chrom_id = chrom_map[chrom_str]
                logger.info(f"Processing {chrom_str}...")

                try:
                    for row in tabix.fetch(chrom_str):
                        try:
                            fields = row.strip().split("\t")
                            if len(fields) < 4:
                                logger.warning(
                                    f"Invalid fragment format: expected at least 4 columns, got {len(fields)}"
                                )
                                continue

                            cell = fields[3]
                            if cell not in valid_barcodes:
                                continue

                            found_cells.add(cell)  # Track that we found this cell

                            start = int(fields[1])
                            end = int(fields[2])

                            if start < 0 or end < 0 or start >= end:
                                logger.warning(f"Invalid fragment coordinates: start={start}, end={end}")
                                continue

                            # Determine bin range spanned by this fragment
                            bin_start = start // self.bin_size
                            bin_end = (end - 1) // self.bin_size  # ensure end is inclusive if needed

                            cell_type = cell_type_map[cell]

                            # Count both Tn5 insertion sites independently for ATAC-seq accessibility
                            # Fragment intervals represent accessible chromatin between insertion sites
                            # See: https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/outputs/fragments-file#fragment-interval-5b7699
                            coverage_aggregator[(chrom_id, bin_start, cell_type)] += 1  # Start insertion site
                            if bin_start != bin_end:  # Only add end bin if in different bin
                                coverage_aggregator[(chrom_id, bin_end, cell_type)] += 1  # End insertion site

                        except ValueError as e:
                            logger.warning(f"Failed to parse fragment row '{row.strip()}': {e}")
                            continue
                except ValueError:
                    continue

        # Log missing cells at the end
        missing_cells = len(valid_barcodes - found_cells)
        if missing_cells:
            logger.warning(f"{missing_cells} barcodes in .h5ad were not found in fragment file.")

        if not coverage_aggregator:
            return

        # Convert aggregated data to DataFrame for final processing
        data_tuples = [
            (chrom, bin_id, cell_type, count) for (chrom, bin_id, cell_type), count in coverage_aggregator.items()
        ]

        full_df = pd.DataFrame(data_tuples, columns=["chrom", "bin", "cell_type", "coverage"])

        # Compute total coverage across all chromosomes
        cell_type_totals = full_df.groupby("cell_type")["coverage"].sum()

        full_df["total_coverage"] = full_df["cell_type"].map(cell_type_totals)
        full_df["normalized_coverage"] = (
            (full_df["coverage"] / full_df["total_coverage"]) * self.normalization_factor
        ).fillna(0)

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

    def process_fragment_file(
        self, obs: pd.DataFrame, array_name: str
    ) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
        df_meta, genome_version = self.extract_cell_metadata_from_h5ad(obs)
        valid_barcodes = set(df_meta["cell_name"])
        cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))

        max_chrom, chrom_map = self.build_chrom_mapping(genome_version)

        max_bins = self.calculate_max_bins(genome_version)
        self.create_dataframe_array(array_name, max_chrom, max_bins)
        self.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

        # Build cell_id_map for return value (maintaining API compatibility)
        cell_id_map = {cell: {"cell_type": cell_type_map[cell]} for cell in valid_barcodes}

        return df_meta, cell_id_map
