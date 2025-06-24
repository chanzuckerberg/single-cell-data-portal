import logging
import os
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd
import pysam
import tiledb
from tqdm import tqdm

from .config.genome_config import CHROM_LENGTHS

logger = logging.getLogger(__name__)


def _process_chromosome_worker(args):
    """Worker function for parallel chromosome processing."""
    fragment_file, chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size = args

    coverage_aggregator = defaultdict(int)
    found_cells = set()

    try:
        with pysam.TabixFile(fragment_file) as tabix:
            for row in tabix.fetch(chrom_str):
                try:
                    fields = row.strip().split("\t")
                    if len(fields) < 4:
                        continue

                    cell = fields[3]
                    if cell not in valid_barcodes:
                        continue

                    found_cells.add(cell)

                    start, end = int(fields[1]), int(fields[2])
                    if start < 0 or end < 0 or start >= end:
                        continue

                    # Calculate bins and update coverage
                    bin_start = start // bin_size
                    bin_end = (end - 1) // bin_size
                    cell_type = cell_type_map[cell]

                    # Count both Tn5 insertion sites independently
                    coverage_aggregator[(chrom_id, bin_start, cell_type)] += 1
                    if bin_start != bin_end:
                        coverage_aggregator[(chrom_id, bin_end, cell_type)] += 1

                except ValueError:
                    continue

    except ValueError:
        # Handle case where chromosome doesn't exist in tabix file
        pass

    return dict(coverage_aggregator), found_cells


# Config
BIN_SIZE = 100
NORMALIZATION_FACTOR = 2_000_000


class ATACDataProcessor:
    def __init__(
        self,
        fragment_artifact_id: Optional[str] = None,
        ctx: Optional[tiledb.Ctx] = None,
        enable_parallel: bool = True,
    ) -> None:
        if fragment_artifact_id is not None and not os.path.exists(fragment_artifact_id):
            raise FileNotFoundError(f"Fragment file not found: {fragment_artifact_id}")

        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.chrom_lengths = CHROM_LENGTHS
        self.normalization_factor = NORMALIZATION_FACTOR
        self.enable_parallel = enable_parallel

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
        """Main orchestration method for processing fragment data and writing to TileDB."""
        # Process all fragments and aggregate coverage
        coverage_aggregator, found_cells = self._process_all_chromosomes(
            chrom_map, cell_type_map, valid_barcodes, self.enable_parallel
        )

        # Report missing cells
        self._report_missing_cells(valid_barcodes, found_cells)

        # Early exit if no coverage data found (empty fragment file, no matching cells, or processing errors)
        if not coverage_aggregator:
            return

        # Convert to DataFrame and normalize
        coverage_df = self._create_coverage_dataframe(coverage_aggregator)

        # Write to TileDB
        self._write_coverage_to_tiledb(array_name, coverage_df)

    def _process_all_chromosomes(
        self,
        chrom_map: Dict[str, int],
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
        enable_parallel: bool = True,
    ) -> Tuple[defaultdict, Set[str]]:
        """Process fragments from all chromosomes and aggregate coverage counts using parallel processing."""
        coverage_aggregator = defaultdict(int)  # (chrom, bin, cell_type) -> coverage_count
        found_cells = set()

        # Determine number of processes to use
        num_processes = min(cpu_count(), len(chrom_map)) if enable_parallel else 1

        if num_processes > 1 and enable_parallel:
            logger.info(f"Using {num_processes} processes for parallel chromosome processing")

            # Prepare arguments for parallel processing
            chrom_args = [
                (
                    self.fragment_artifact_id,
                    chrom_str,
                    chrom_map[chrom_str],
                    cell_type_map,
                    valid_barcodes,
                    self.bin_size,
                )
                for chrom_str in chrom_map
            ]

            # Process chromosomes in parallel
            with Pool(processes=num_processes) as pool:
                results = list(
                    tqdm(
                        pool.imap(_process_chromosome_worker, chrom_args),
                        total=len(chrom_args),
                        desc="Processing chromosomes",
                        unit="chrom",
                    )
                )

            # Aggregate results from all processes
            for chrom_coverage, chrom_found_cells in results:
                for key, count in chrom_coverage.items():
                    coverage_aggregator[key] += count
                found_cells.update(chrom_found_cells)
        else:
            # Fallback to serial processing
            logger.info("Using serial processing (single process)")
            with pysam.TabixFile(self.fragment_artifact_id) as tabix:
                for chrom_str in tqdm(chrom_map, desc="Processing chromosomes", unit="chrom"):
                    chrom_id = chrom_map[chrom_str]
                    logger.info(f"Processing {chrom_str}...")
                    self._process_chromosome(
                        tabix, chrom_str, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
                    )

        return coverage_aggregator, found_cells

    def _process_chromosome(
        self,
        tabix: pysam.TabixFile,
        chrom_str: str,
        chrom_id: int,
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
        coverage_aggregator: defaultdict,
        found_cells: Set[str],
    ) -> None:
        """Process fragments for a single chromosome."""
        try:
            for row in tabix.fetch(chrom_str):
                self._process_fragment_row(
                    row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
                )
        except ValueError:
            logger.warning(f"Failed to fetch chromosome {chrom_str}")

    def _process_fragment_row(
        self,
        row: str,
        chrom_id: int,
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
        coverage_aggregator: defaultdict,
        found_cells: Set[str],
    ) -> None:
        """Process a single fragment row and update coverage counts."""
        try:
            fields = row.strip().split("\t")
            if len(fields) < 4:
                logger.warning(f"Invalid fragment format: expected at least 4 columns, got {len(fields)}")
                return

            cell = fields[3]
            if cell not in valid_barcodes:
                return

            found_cells.add(cell)

            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end < 0 or start >= end:
                logger.warning(f"Invalid fragment coordinates: start={start}, end={end}")
                return

            # Calculate bins and update coverage
            bin_start = start // self.bin_size
            bin_end = (end - 1) // self.bin_size
            cell_type = cell_type_map[cell]

            # Count both Tn5 insertion sites independently for ATAC-seq accessibility
            # Fragment intervals represent accessible chromatin between insertion sites
            # See: https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/outputs/fragments-file#fragment-interval-5b7699
            coverage_aggregator[(chrom_id, bin_start, cell_type)] += 1  # Start insertion site
            if bin_start != bin_end:  # Only add end bin if in different bin
                coverage_aggregator[(chrom_id, bin_end, cell_type)] += 1  # End insertion site

        except ValueError as e:
            logger.warning(f"Failed to parse fragment row '{row.strip()}': {e}")

    def _report_missing_cells(self, valid_barcodes: Set[str], found_cells: Set[str]) -> None:
        """Log warning about cells not found in fragment file."""
        missing_cells = len(valid_barcodes - found_cells)
        if missing_cells:
            logger.warning(f"{missing_cells} barcodes in .h5ad were not found in fragment file.")

    def _compute_cell_type_totals(self, coverage_aggregator: defaultdict) -> Dict[str, int]:
        """Compute total coverage per cell type for normalization."""
        logger.info("Computing cell type totals for normalization...")
        cell_type_totals = defaultdict(int)
        for (_, _, cell_type), count in coverage_aggregator.items():
            cell_type_totals[cell_type] += count
        return dict(cell_type_totals)

    def _process_coverage_data(
        self,
        coverage_aggregator: defaultdict,
    ) -> None:
        """Process coverage data for DataFrame creation."""
        cell_type_totals = self._compute_cell_type_totals(coverage_aggregator)
        total_records = len(coverage_aggregator)

        logger.info(f"Processing {total_records:,} records for DataFrame creation...")

        for _, ((chrom, bin_id, cell_type), count) in enumerate(
            tqdm(coverage_aggregator.items(), desc="Processing DataFrame creation", unit="records")
        ):
            total_coverage = cell_type_totals[cell_type]
            normalized_coverage = (count / total_coverage) * self.normalization_factor if total_coverage > 0 else 0.0

            record = {
                "chrom": chrom,
                "bin": bin_id,
                "cell_type": cell_type,
                "coverage": count,
                "total_coverage": total_coverage,
                "normalized_coverage": normalized_coverage,
            }

            # Add record to current chunk
            self._current_chunk.append(record)

            # Process chunk when it reaches target size
            if len(self._current_chunk) >= self._dataframe_chunk_size:
                self._chunks.append(pd.DataFrame(self._current_chunk))
                self._current_chunk.clear()

    def _create_coverage_dataframe(self, coverage_aggregator: defaultdict, chunk_size: int = 50000) -> pd.DataFrame:
        """Convert aggregated coverage data to normalized DataFrame using chunked processing."""
        # Store as instance variables temporarily
        self._chunks = []
        self._current_chunk = []
        self._dataframe_chunk_size = chunk_size

        # Process coverage data for DataFrame creation
        self._process_coverage_data(coverage_aggregator)

        # Process remaining records
        if self._current_chunk:
            self._chunks.append(pd.DataFrame(self._current_chunk))

        logger.info("Concatenating chunks...")
        if not self._chunks:
            # Clean up instance variables
            del self._chunks
            del self._current_chunk
            del self._dataframe_chunk_size
            return pd.DataFrame(
                columns=["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
            )

        # Concatenate all chunks efficiently
        df = pd.concat(self._chunks, ignore_index=True, copy=False)

        # Clean up instance variables
        del self._chunks
        del self._current_chunk
        del self._dataframe_chunk_size

        # Optimize data types to save memory
        df["chrom"] = df["chrom"].astype("int32")
        df["bin"] = df["bin"].astype("int32")
        df["coverage"] = df["coverage"].astype("int32")
        df["total_coverage"] = df["total_coverage"].astype("int32")
        df["normalized_coverage"] = df["normalized_coverage"].astype("float32")

        logger.info(f"Created DataFrame with {len(df):,} records")
        return df

    def _write_coverage_to_tiledb(self, array_name: str, coverage_df: pd.DataFrame) -> None:
        """Write coverage DataFrame to TileDB array."""
        logger.info(f"Writing {len(coverage_df):,} coverage records to TileDB...")
        with tiledb.SparseArray(array_name, mode="w", ctx=self.ctx) as A:
            A[
                (
                    coverage_df["chrom"].astype("int32").to_numpy(),
                    coverage_df["bin"].astype("int32").to_numpy(),
                    coverage_df["cell_type"].to_numpy(),
                )
            ] = {
                "coverage": coverage_df["coverage"].astype("int32").to_numpy(),
                "total_coverage": coverage_df["total_coverage"].astype("int32").to_numpy(),
                "normalized_coverage": coverage_df["normalized_coverage"].to_numpy(),
            }
        logger.info("Successfully wrote coverage data to TileDB")

    def process_fragment_file(self, obs: pd.DataFrame, array_name: str) -> None:
        df_meta, genome_version = self.extract_cell_metadata_from_h5ad(obs)

        valid_barcodes = set(df_meta["cell_name"])
        cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))

        max_chrom, chrom_map = self.build_chrom_mapping(genome_version)
        max_bins = self.calculate_max_bins(genome_version)

        self.create_dataframe_array(array_name, max_chrom, max_bins)
        self.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)
