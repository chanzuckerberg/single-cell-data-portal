import itertools
import logging
import os
import tempfile
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Optional, Set, Tuple

import boto3
import cellxgene_schema.atac_seq as atac_seq
import numpy as np
import pandas as pd
import pysam
import tiledb
from tqdm import tqdm

logger = logging.getLogger(__name__)


BIN_SIZE = 100
NORMALIZATION_FACTOR = 2_000_000


class ATACDataProcessor:
    def __init__(
        self,
        fragment_artifact_id: Optional[str] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ) -> None:
        if fragment_artifact_id is not None and not self._file_exists(fragment_artifact_id):
            raise FileNotFoundError(f"Fragment file not found: {fragment_artifact_id}")

        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.normalization_factor = NORMALIZATION_FACTOR
        self._local_fragment_file = None
        self._local_fragment_index = None

    def _file_exists(self, file_path: str) -> bool:
        """Check if a file exists, supporting both local paths and S3 URLs."""
        if file_path.startswith("s3://"):
            s3_path = file_path[5:]
            bucket, key = s3_path.split("/", 1)

            try:
                s3_client = boto3.client("s3")
                s3_client.head_object(Bucket=bucket, Key=key)
                return True
            except Exception:
                return False
        else:
            return os.path.exists(file_path)

    def _download_s3_file(self, s3_path: str) -> str:
        """Download S3 file and its index to local temporary files and return local path."""
        if not s3_path.startswith("s3://"):
            return s3_path

        s3_url = s3_path[5:]
        bucket, key = s3_url.split("/", 1)

        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".bgz")
        temp_path = temp_file.name
        temp_file.close()

        temp_index_path = temp_path + ".tbi"

        try:
            logger.info(f"Downloading S3 file {s3_path} to {temp_path}")
            s3_client = boto3.client("s3")

            s3_client.download_file(bucket, key, temp_path)
            logger.info(f"Successfully downloaded S3 file to {temp_path}")

            index_key = key + ".tbi"
            logger.info(f"Downloading S3 index file s3://{bucket}/{index_key} to {temp_index_path}")
            s3_client.download_file(bucket, index_key, temp_index_path)
            logger.info(f"Successfully downloaded S3 index file to {temp_index_path}")

            return temp_path
        except Exception as e:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            if os.path.exists(temp_index_path):
                os.unlink(temp_index_path)
            raise RuntimeError(f"Failed to download S3 file {s3_path}: {str(e)}") from e

    def _get_fragment_file_path(self) -> str:
        """Get the fragment file path, downloading from S3 if necessary."""
        if self.fragment_artifact_id is None:
            raise ValueError("Fragment artifact ID is not set")

        if self._local_fragment_file is None:
            self._local_fragment_file = self._download_s3_file(self.fragment_artifact_id)
            if self._local_fragment_file != self.fragment_artifact_id:
                self._local_fragment_index = self._local_fragment_file + ".tbi"

        return self._local_fragment_file

    def _cleanup_temp_files(self) -> None:
        """Clean up temporary downloaded files."""
        if (
            self._local_fragment_file
            and self._local_fragment_file != self.fragment_artifact_id
            and os.path.exists(self._local_fragment_file)
        ):
            logger.info(f"Cleaning up temporary file: {self._local_fragment_file}")
            os.unlink(self._local_fragment_file)
            self._local_fragment_file = None

        if self._local_fragment_index and os.path.exists(self._local_fragment_index):
            logger.info(f"Cleaning up temporary index file: {self._local_fragment_index}")
            os.unlink(self._local_fragment_index)
            self._local_fragment_index = None

    @staticmethod
    def _process_chromosome_worker(args):
        """Worker function for parallel chromosome processing."""
        fragment_file, chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size = args

        coverage_aggregator = defaultdict(int)
        found_cells = set()

        try:
            with pysam.TabixFile(fragment_file) as tabix:
                fragment_iter = tabix.fetch(chrom_str)

                for row in fragment_iter:
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

                        bin_start = start // bin_size
                        bin_end = (end - 1) // bin_size
                        cell_type = cell_type_map[cell]

                        coverage_aggregator[(chrom_id, bin_start, cell_type)] += 1
                        if bin_start != bin_end:
                            coverage_aggregator[(chrom_id, bin_end, cell_type)] += 1

                    except ValueError:
                        continue

        except ValueError:
            pass

        return dict(coverage_aggregator), found_cells

    def extract_cell_metadata_from_h5ad(
        self, obs: pd.DataFrame, obs_column: str = "cell_type", uns: Optional[Dict] = None
    ) -> Tuple[pd.DataFrame, str]:
        if obs_column not in obs.columns:
            raise ValueError(f"Column {obs_column} not found in obs DataFrame.")

        organism_ontology_term_id = None

        if uns is not None and "organism_ontology_term_id" in uns:
            organism_ontology_term_id = uns["organism_ontology_term_id"]
        elif "organism_ontology_term_id" in obs.columns:
            organism_ontology_term_id = obs["organism_ontology_term_id"].iloc[0]
        else:
            raise ValueError(
                "Column 'organism_ontology_term_id' is required but not found in obs DataFrame or uns dict."
            )

        if pd.isna(organism_ontology_term_id):
            raise ValueError("organism_ontology_term_id cannot be null/NaN.")

        df = obs[[obs_column]].copy()
        df = df.rename_axis("cell_name").reset_index()
        chromosome_by_length = atac_seq.organism_ontology_term_id_by_chromosome_length_table.get(
            organism_ontology_term_id
        )
        return df, chromosome_by_length

    def calculate_max_bins(self, chromosome_by_length: Dict) -> int:
        bins_per_chrom = {chrom: int(length / self.bin_size) + 1 for chrom, length in chromosome_by_length.items()}
        return max(bins_per_chrom.values()) - 1

    def build_chrom_mapping(self, chromosome_by_length: Dict) -> Tuple[int, Dict[str, int]]:
        chrom_map = defaultdict(lambda: 0)
        for i, chrom in enumerate(chromosome_by_length, 1):
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
        self,
        array_name: str,
        chrom_map: Dict[str, int],
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
    ) -> None:
        """Main orchestration method for processing fragment data and writing to TileDB."""
        num_chromosomes = len(chrom_map)
        chromosome_batch_size = cpu_count()

        if num_chromosomes > chromosome_batch_size:
            logger.info(
                f"Using batch processing for {num_chromosomes} chromosomes with batch size {chromosome_batch_size}"
            )
            self._process_chromosomes_in_batches(
                array_name, chrom_map, cell_type_map, valid_barcodes, chromosome_batch_size
            )
        else:
            logger.info(f"Using single-pass processing for {num_chromosomes} chromosomes")
            coverage_aggregator, found_cells = self._process_all_chromosomes(chrom_map, cell_type_map, valid_barcodes)

            self._report_missing_cells(valid_barcodes, found_cells)

            if not coverage_aggregator:
                return

            # For single-pass processing, compute totals directly
            cell_type_totals = self._compute_cell_type_totals_from_aggregator(coverage_aggregator)
            self._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)

    def _process_all_chromosomes(
        self,
        chrom_map: Dict[str, int],
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
    ) -> Tuple[defaultdict, Set[str]]:
        """Process fragments from all chromosomes and aggregate coverage counts."""
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        num_processes = min(cpu_count(), max(1, len(chrom_map)))
        logger.info(f"Using {num_processes} processes for chromosome processing")

        fragment_file_path = self._get_fragment_file_path()

        chrom_args = [
            (
                fragment_file_path,
                chrom_str,
                chrom_map[chrom_str],
                cell_type_map,
                valid_barcodes,
                self.bin_size,
            )
            for chrom_str in chrom_map
        ]

        with Pool(processes=num_processes) as pool:
            results = list(
                tqdm(
                    pool.imap(self._process_chromosome_worker, chrom_args),
                    total=len(chrom_args),
                    desc="Processing chromosomes",
                    unit="chrom",
                )
            )

        for chrom_coverage, chrom_found_cells in results:
            for key, count in chrom_coverage.items():
                coverage_aggregator[key] += count
            found_cells.update(chrom_found_cells)

        return coverage_aggregator, found_cells

    def _report_missing_cells(self, valid_barcodes: Set[str], found_cells: Set[str]) -> None:
        """Log warning about cells not found in fragment file."""
        missing_cells = len(valid_barcodes - found_cells)
        if missing_cells:
            logger.warning(f"{missing_cells} barcodes in .h5ad were not found in fragment file.")

    def _process_chromosomes_in_batches(
        self,
        array_name: str,
        chrom_map: Dict[str, int],
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
        chromosome_batch_size: int,
    ) -> None:
        """Process chromosomes in batches to reduce memory usage."""
        all_found_cells = set()
        all_coverage_aggregator = defaultdict(int)

        chrom_list = list(chrom_map.items())
        chrom_batches = []
        for i in range(0, len(chrom_list), chromosome_batch_size):
            batch = chrom_list[i : i + chromosome_batch_size]
            chrom_batches.append(batch)

        logger.info(f"Processing {len(chrom_list)} chromosomes in {len(chrom_batches)} batches")

        logger.info("Computing global cell type totals for normalization...")
        global_cell_type_totals = self._compute_global_cell_type_totals(chrom_batches, cell_type_map, valid_barcodes)

        for batch_idx, chrom_batch in enumerate(chrom_batches):
            logger.info(f"Processing batch {batch_idx + 1}/{len(chrom_batches)} with {len(chrom_batch)} chromosomes")

            batch_chrom_map = dict(chrom_batch)
            batch_coverage_aggregator, batch_found_cells = self._process_all_chromosomes(
                batch_chrom_map, cell_type_map, valid_barcodes
            )

            all_found_cells.update(batch_found_cells)

            if batch_coverage_aggregator:
                for key, count in batch_coverage_aggregator.items():
                    all_coverage_aggregator[key] += count
            else:
                logger.info(f"Batch {batch_idx + 1} had no coverage data, skipping")

        self._report_missing_cells(valid_barcodes, all_found_cells)

        if all_coverage_aggregator:
            self._stream_coverage_chunks_to_tiledb(all_coverage_aggregator, global_cell_type_totals, array_name)

    def _compute_global_cell_type_totals(
        self,
        chrom_batches: List[List[Tuple[str, int]]],
        cell_type_map: Dict[str, str],
        valid_barcodes: Set[str],
    ) -> Dict[str, int]:
        """Compute global cell type totals across all batches for normalization."""
        global_cell_type_totals = defaultdict(int)

        for batch_idx, chrom_batch in enumerate(chrom_batches):
            logger.info(f"Computing totals for batch {batch_idx + 1}/{len(chrom_batches)}")

            batch_chrom_map = dict(chrom_batch)
            batch_coverage_aggregator, _ = self._process_all_chromosomes(batch_chrom_map, cell_type_map, valid_barcodes)

            for (_, _, cell_type), count in batch_coverage_aggregator.items():
                global_cell_type_totals[cell_type] += count

        return dict(global_cell_type_totals)

    def _compute_cell_type_totals_from_aggregator(self, coverage_aggregator: defaultdict) -> Dict[str, int]:
        """Compute cell type totals from coverage aggregator data."""
        cell_type_totals = defaultdict(int)
        for (_, _, cell_type), count in coverage_aggregator.items():
            cell_type_totals[cell_type] += count
        return dict(cell_type_totals)

    def _stream_coverage_chunks_to_tiledb(
        self,
        coverage_aggregator: defaultdict,
        global_cell_type_totals: Dict[str, int],
        array_name: str,
        chunk_size: int = 50000,
    ) -> int:
        """Streaming approach: process directly into pre-allocated numpy arrays, write once."""

        total_records = len(coverage_aggregator)
        logger.info(f"Processing {total_records:,} records with direct array assignment...")
        array_index = 0

        def generate_coverage_records():
            """Generator to yield coverage records for direct array assignment."""
            for (chrom, bin_id, cell_type), count in coverage_aggregator.items():
                total_coverage = global_cell_type_totals.get(cell_type, 0)
                normalized_coverage = (
                    (count / total_coverage) * self.normalization_factor if total_coverage > 0 else 0.0
                )
                yield chrom, bin_id, cell_type, count, total_coverage, normalized_coverage

        coverage_records = generate_coverage_records()
        while True:
            chunk = list(itertools.islice(coverage_records, chunk_size))
            array_index += len(chunk)
            if not chunk:
                break
            write_range = f"{array_index-len(chunk):,}-{array_index:,}"
            logger.info(f"Writing {write_range} records to TileDB as single fragment...")
            self._write_arrays_to_tiledb(array_name, *zip(*chunk, strict=True))
        return array_index

    def _write_arrays_to_tiledb(
        self,
        array_name: str,
        chroms: np.ndarray,
        bins: np.ndarray,
        cell_types: np.ndarray,
        coverages: np.ndarray,
        total_coverages: np.ndarray,
        normalized_coverages: np.ndarray,
    ) -> None:
        """Write numpy arrays directly to TileDB for maximum efficiency."""
        logger.info(f"Writing {len(chroms):,} records to TileDB from numpy arrays...")
        with tiledb.SparseArray(array_name, mode="w", ctx=self.ctx) as A:
            A[(chroms, bins, cell_types)] = {
                "coverage": coverages,
                "total_coverage": total_coverages,
                "normalized_coverage": normalized_coverages,
            }
        logger.info("Successfully wrote numpy arrays to TileDB")

    def process_fragment_file(self, obs: pd.DataFrame, array_name: str, uns: Optional[Dict] = None) -> None:
        try:
            logger.info(f"Starting ATAC fragment processing for array: {array_name}")
            df_meta, chromosome_by_length = self.extract_cell_metadata_from_h5ad(obs, uns=uns)

            valid_barcodes = set(df_meta["cell_name"])
            cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))

            max_chrom, chrom_map = self.build_chrom_mapping(chromosome_by_length)
            max_bins = self.calculate_max_bins(chromosome_by_length)

            logger.info(f"Creating TileDB array with {max_chrom} chromosomes and {max_bins} max bins")
            self.create_dataframe_array(array_name, max_chrom, max_bins)

            logger.info(f"Processing {len(chrom_map)} chromosomes with {len(valid_barcodes)} valid barcodes")
            self.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

            logger.info("ATAC fragment processing completed successfully")

        except Exception as e:
            logger.error(f"ATAC fragment processing failed: {e}")
            raise

        finally:
            self._cleanup_temp_files()
