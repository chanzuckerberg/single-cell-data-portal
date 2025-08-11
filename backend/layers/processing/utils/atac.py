import itertools
import logging
import os
import tempfile
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from typing import Dict, Optional, Set, Tuple

import boto3
import cellxgene_schema.atac_seq as atac_seq
import numpy as np
import pandas as pd
import psutil
import pysam
import tiledb
from tqdm import tqdm

logger = logging.getLogger(__name__)


BIN_SIZE = 100
NORMALIZATION_FACTOR = 2_000_000


def get_dynamic_chunk_size(memory_factor: float = 0.1, min_chunk_size: int = 100_000, max_chunk_size: int = 50_000_000):
    """
    Calculate dynamic chunk size based on available memory.

    Args:
        memory_factor: Fraction of available memory to use for chunk sizing
        min_chunk_size: Minimum chunk size regardless of memory
        max_chunk_size: Maximum chunk size to prevent excessive memory usage

    Returns:
        Optimal chunk size for the current system
    """
    available_memory = psutil.virtual_memory().available

    # Use available memory to be more conservative
    memory_for_chunks = int(available_memory * memory_factor)

    # Estimate memory per record (for ATAC data):
    # - chromosome name: ~10 bytes
    # - bin_id: 8 bytes (int64)
    # - cell_type: ~20 bytes average
    # - coverage: 8 bytes (int64)
    # - total_coverage: 8 bytes (int64)
    # - normalized_coverage: 8 bytes (float64)
    # - dict overhead: ~100 bytes
    # Total: ~162 bytes per record
    bytes_per_record = 162

    calculated_chunk_size = memory_for_chunks // bytes_per_record

    # Apply min/max constraints
    chunk_size = max(min_chunk_size, min(calculated_chunk_size, max_chunk_size))

    return chunk_size


class ATACDataProcessor:
    def __init__(
        self,
        fragment_artifact_id: Optional[str] = None,
        ctx: Optional[tiledb.Ctx] = None,
        min_coverage_threshold: float = 0.1,  # Set to 0 to disable pruning entirely
        enable_quantization: bool = True,  # Enable normalized coverage quantization for better compression
    ) -> None:
        if fragment_artifact_id is not None and not self._file_exists(fragment_artifact_id):
            raise FileNotFoundError(f"Fragment file not found: {fragment_artifact_id}")

        self.fragment_artifact_id = fragment_artifact_id
        self.ctx = ctx
        self.bin_size = BIN_SIZE
        self.normalization_factor = NORMALIZATION_FACTOR
        self.min_coverage_threshold = min_coverage_threshold
        self.enable_quantization = enable_quantization
        self._local_fragment_file = None
        self._local_fragment_index = None

    def _quantize_normalized_coverage(self, normalized_coverage: float) -> int:
        """
        Quantize normalized coverage to uint8 for better compression.

        Uses linear quantization mapping 0-10 range to 0-255 (uint8).
        Precision: ~0.04 units, Mean error: ~0.02 units.
        Provides ~6x total compression improvement.

        Args:
            normalized_coverage: Float value to quantize

        Returns:
            Quantized integer value (0-255)
        """
        if not self.enable_quantization:
            return normalized_coverage

        # Linear quantization: 0-10 range -> 0-255 (uint8)
        # Clip extreme values to prevent overflow
        quantized = np.clip(normalized_coverage * 25.5, 0, 255)
        return int(quantized)

    def _dequantize_normalized_coverage(self, quantized_value: int) -> float:
        """
        Convert quantized uint8 back to normalized coverage.

        Args:
            quantized_value: Quantized integer (0-255)

        Returns:
            Reconstructed float value
        """
        if not self.enable_quantization:
            return float(quantized_value)

        # Reverse linear quantization: 0-255 -> 0-10 range
        return quantized_value / 25.5

    def _file_exists(self, file_path: str) -> bool:
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
        coverage_int_compression = tiledb.FilterList(
            [
                tiledb.BitWidthReductionFilter(),
                tiledb.ByteShuffleFilter(),
                tiledb.ZstdFilter(level=6),
            ]
        )

        coverage_float_compression = tiledb.FilterList(
            [
                tiledb.ByteShuffleFilter(),
                tiledb.ZstdFilter(level=6),
            ]
        )

        categorical_compression = tiledb.FilterList(
            [
                tiledb.DictionaryFilter(),
                tiledb.ZstdFilter(level=22),
            ]
        )

        dim_filters = tiledb.FilterList(
            [
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(level=6),
            ]
        )

        # Adaptive tile sizing: optimize for dataset size and query patterns
        calculated_tile = min(max(max_bins // 1000, 100), 10000)
        optimal_bin_tile = min(calculated_tile, max_bins)
        genomic_window_kb = (optimal_bin_tile * self.bin_size) // 1000

        logger.info(
            f"Adaptive tiling: {optimal_bin_tile} bins per tile ({genomic_window_kb}kb genomic windows) "
            f"for dataset with {max_bins:,} total bins"
        )

        domain = tiledb.Domain(
            tiledb.Dim(name="chrom", domain=(1, max_chrom), tile=1, dtype=np.uint32, filters=dim_filters),
            tiledb.Dim(name="bin", domain=(0, max_bins), tile=optimal_bin_tile, dtype=np.uint32, filters=dim_filters),
            tiledb.Dim(name="cell_type", dtype="ascii", filters=categorical_compression),
        )

        logger.info(
            "Using advanced compression filters: BitWidthReduction+ByteShuffle+ZStd for integers, "
            "DoubleDelta+BitWidthReduction+ZStd for coordinates"
        )

        # Choose storage type based on quantization setting
        if self.enable_quantization:
            normalized_dtype = np.uint8
            normalized_filters = coverage_int_compression
            logger.info("Using quantized normalized coverage (uint8) for ~6x better compression")
        else:
            normalized_dtype = np.float32
            normalized_filters = coverage_float_compression

        schema = tiledb.ArraySchema(
            domain=domain,
            attrs=[
                tiledb.Attr(name="coverage", dtype=np.uint16, filters=coverage_int_compression),
                tiledb.Attr(name="total_coverage", dtype=np.uint32, filters=coverage_int_compression),
                tiledb.Attr(name="normalized_coverage", dtype=normalized_dtype, filters=normalized_filters),
            ],
            sparse=True,
            allows_duplicates=False,
        )
        tiledb.SparseArray.create(array_name, schema)
        logger.info(
            f"Created TileDB array {array_name} with advanced filter compression: "
            f"BitWidth+ByteShuffle+ZStd(level=6) for coverage, DoubleDelta+BitWidth+ZStd(level=6) for coordinates, "
            f"Dictionary+ZStd(level=22) for cell types"
        )

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

        # Memory-aware process limiting to prevent OOM
        fragment_file_path = self._get_fragment_file_path()
        fragment_file_size_gb = os.path.getsize(fragment_file_path) / (1024**3)
        available_memory_gb = psutil.virtual_memory().total * 0.7 / (1024**3)  # Use 70% of available memory
        max_memory_processes = max(1, int(available_memory_gb / fragment_file_size_gb))

        num_processes = min(cpu_count(), len(chrom_map), max_memory_processes)
        logger.info(
            f"Fragment file size: {fragment_file_size_gb:.2f}GB, using {num_processes} processes (memory-limited from {cpu_count()})"
        )

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

        logger.info("Processing chromosomes with single-pass totals computation...")
        global_cell_type_totals = defaultdict(int)

        for batch_idx, chrom_batch in enumerate(chrom_batches):
            logger.info(f"Processing batch {batch_idx + 1}/{len(chrom_batches)} with {len(chrom_batch)} chromosomes")

            batch_chrom_map = dict(chrom_batch)
            batch_coverage_aggregator, batch_found_cells = self._process_all_chromosomes(
                batch_chrom_map, cell_type_map, valid_barcodes
            )

            all_found_cells.update(batch_found_cells)

            if batch_coverage_aggregator:
                # Compute totals incrementally to avoid double-pass
                for key, count in batch_coverage_aggregator.items():
                    _, _, cell_type = key
                    global_cell_type_totals[cell_type] += count
                    all_coverage_aggregator[key] += count
            else:
                logger.info(f"Batch {batch_idx + 1} had no coverage data, skipping")

        self._report_missing_cells(valid_barcodes, all_found_cells)

        if all_coverage_aggregator:
            self._stream_coverage_chunks_to_tiledb(all_coverage_aggregator, dict(global_cell_type_totals), array_name)

    def _compute_cell_type_totals_from_aggregator(self, coverage_aggregator: defaultdict) -> Dict[str, int]:
        cell_type_totals = defaultdict(int)
        for (_, _, cell_type), count in coverage_aggregator.items():
            cell_type_totals[cell_type] += count
        return dict(cell_type_totals)

    def _stream_coverage_chunks_to_tiledb(
        self,
        coverage_aggregator: defaultdict,
        global_cell_type_totals: Dict[str, int],
        array_name: str,
        chunk_size: int = None,
    ) -> int:
        """Generator approach: process chunks on-the-fly within single TileDB session."""

        if chunk_size is None:
            chunk_size = get_dynamic_chunk_size()

        total_records = len(coverage_aggregator)
        logger.info(
            f"Processing {total_records:,} records with generator-based chunked streaming (chunk_size: {chunk_size:,})..."
        )

        records_processed = self._write_chunks_generator_to_tiledb(
            array_name, self._generate_chunks(coverage_aggregator, global_cell_type_totals, chunk_size)
        )

        logger.info(f"Successfully processed {records_processed:,} records to TileDB as single fragment")
        return records_processed

    def _generate_chunks(
        self, coverage_aggregator: defaultdict, global_cell_type_totals: Dict[str, int], chunk_size: int
    ):
        total_records = len(coverage_aggregator)
        pruned_records = 0
        filtered_items = []

        for (chrom, bin_id, cell_type), count in coverage_aggregator.items():
            total_coverage = global_cell_type_totals.get(cell_type, 0)
            normalized_coverage = (count / total_coverage) * self.normalization_factor if total_coverage > 0 else 0.0

            # Apply sparse data pruning - skip bins with normalized coverage below threshold
            if normalized_coverage < self.min_coverage_threshold:
                pruned_records += 1
                continue

            if self.enable_quantization:
                stored_normalized_coverage = self._quantize_normalized_coverage(normalized_coverage)
            else:
                stored_normalized_coverage = normalized_coverage

            filtered_items.append(
                {
                    "chrom": chrom,
                    "bin_id": bin_id,
                    "cell_type": cell_type,
                    "coverage": count,
                    "total_coverage": total_coverage,
                    "normalized_coverage": stored_normalized_coverage,
                }
            )

        # Log pruning statistics
        if total_records > 0:
            pruning_percent = (pruned_records / total_records) * 100
            kept_records = total_records - pruned_records
            logger.info(
                f"Sparse data pruning (normalized coverage threshold ≥{self.min_coverage_threshold}): "
                f"{pruned_records:,}/{total_records:,} records ({pruning_percent:.1f}%) removed, "
                f"{kept_records:,} records kept"
            )

        items = iter(filtered_items)
        while True:
            chunk = list(itertools.islice(items, chunk_size))
            if not chunk:
                break
            yield chunk

    def _write_chunks_generator_to_tiledb(
        self,
        array_name: str,
        chunk_generator,
    ) -> int:
        records_written = 0
        with tiledb.SparseArray(array_name, mode="w", ctx=self.ctx) as A:
            for chunk_idx, chunk_data in enumerate(tqdm(chunk_generator, desc="Writing chunks to TileDB")):
                chunk_size = len(chunk_data)
                logger.debug(f"Writing chunk {chunk_idx + 1} ({chunk_size:,} records)...")

                chroms = np.array([record["chrom"] for record in chunk_data], dtype=np.int32)
                bins = np.array([record["bin_id"] for record in chunk_data], dtype=np.int32)
                cell_types = np.array([record["cell_type"] for record in chunk_data], dtype=object)

                # Apply optimized data types with overflow protection
                coverage_values = [record["coverage"] for record in chunk_data]
                total_coverage_values = [record["total_coverage"] for record in chunk_data]
                normalized_coverage_values = [record["normalized_coverage"] for record in chunk_data]

                # Check for potential overflows and log warnings
                max_coverage = max(coverage_values) if coverage_values else 0
                max_total_coverage = max(total_coverage_values) if total_coverage_values else 0

                if max_coverage > 65535:  # uint16 max
                    logger.warning(
                        f"Coverage value {max_coverage} exceeds uint16 range (65535), clipping to prevent overflow"
                    )
                    coverage_values = [min(val, 65535) for val in coverage_values]

                if max_total_coverage > 4294967295:  # uint32 max
                    logger.warning(
                        f"Total coverage value {max_total_coverage} exceeds uint32 range (4.3B), clipping to prevent overflow"
                    )
                    total_coverage_values = [min(val, 4294967295) for val in total_coverage_values]

                coverages = np.array(coverage_values, dtype=np.uint16)
                total_coverages = np.array(total_coverage_values, dtype=np.uint32)

                if self.enable_quantization:
                    normalized_coverages = np.array(normalized_coverage_values, dtype=np.uint8)
                else:
                    normalized_coverages = np.array(normalized_coverage_values, dtype=np.float32)

                A[(chroms, bins, cell_types)] = {
                    "coverage": coverages,
                    "total_coverage": total_coverages,
                    "normalized_coverage": normalized_coverages,
                }

                records_written += chunk_size
                logger.debug(f"Successfully wrote chunk {chunk_idx + 1} of {chunk_size:,} records")

        return records_written

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
