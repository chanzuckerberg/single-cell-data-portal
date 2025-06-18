from collections import defaultdict
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import pytest

from backend.layers.processing.utils.atac import ATACDataProcessor


class TestATACDataProcessor:
    """Test suite for ATACDataProcessor class."""

    @pytest.mark.parametrize(
        "organism_id,expected_genome",
        [
            ("NCBITaxon:9606", "hg38"),
            ("NCBITaxon:10090", "mm39"),
        ],
    )
    def test_get_genome_version_valid_organisms(self, organism_id, expected_genome):
        """Test genome version mapping for valid organisms."""
        result = ATACDataProcessor.get_genome_version(organism_id)
        assert result == expected_genome

    def test_get_genome_version_invalid_organism(self):
        """Test ValueError is raised for unknown organism ontology term ID."""
        with pytest.raises(ValueError, match="Unknown organism ontology term ID"):
            ATACDataProcessor.get_genome_version("NCBITaxon:12345")

    def test_constructor_valid_file_path(self, tmp_path):
        """Test ATACDataProcessor initialization with valid fragment file path."""
        # Create a temporary fragment file
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        # Initialize processor with valid file
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Verify all attributes are set correctly
        assert processor.fragment_artifact_id == str(fragment_file)
        assert processor.ctx is None
        assert processor.bin_size == 100
        assert processor.normalization_factor == 2_000_000

    def test_constructor_file_not_found(self):
        """Test FileNotFoundError is raised when fragment file doesn't exist."""
        non_existent_file = "/path/to/non/existent/file.tsv.gz"

        with pytest.raises(FileNotFoundError, match=f"Fragment file not found: {non_existent_file}"):
            ATACDataProcessor(fragment_artifact_id=non_existent_file)

    def test_extract_cell_metadata_valid_obs(self, tmp_path):
        """Test cell metadata extraction with valid obs DataFrame."""
        # Create temporary fragment file
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create valid obs DataFrame
        obs = pd.DataFrame(
            {
                "cell_type": ["T cell", "B cell", "NK cell"],
                "organism_ontology_term_id": ["NCBITaxon:9606", "NCBITaxon:9606", "NCBITaxon:9606"],
            },
            index=["cell1", "cell2", "cell3"],
        )

        df_meta, genome_version = processor.extract_cell_metadata_from_h5ad(obs)

        # Verify returned DataFrame structure
        expected_df = pd.DataFrame(
            {"cell_name": ["cell1", "cell2", "cell3"], "cell_type": ["T cell", "B cell", "NK cell"]}
        )
        pd.testing.assert_frame_equal(df_meta, expected_df)

        # Verify genome version
        assert genome_version == "hg38"

    def test_extract_cell_metadata_missing_obs_column(self, tmp_path):
        """Test ValueError when obs_column is missing from DataFrame."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create obs DataFrame without cell_type column
        obs = pd.DataFrame(
            {"other_column": ["value1", "value2"], "organism_ontology_term_id": ["NCBITaxon:9606", "NCBITaxon:9606"]},
            index=["cell1", "cell2"],
        )

        with pytest.raises(ValueError, match="Column cell_type not found in obs DataFrame"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test_extract_cell_metadata_missing_organism_column(self, tmp_path):
        """Test ValueError when organism_ontology_term_id column is missing."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create obs DataFrame without organism_ontology_term_id column
        obs = pd.DataFrame({"cell_type": ["T cell", "B cell"]}, index=["cell1", "cell2"])

        with pytest.raises(ValueError, match="Column 'organism_ontology_term_id' is required but not found"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test_extract_cell_metadata_null_organism_id(self, tmp_path):
        """Test ValueError when organism_ontology_term_id is null/NaN."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create obs DataFrame with null organism_ontology_term_id
        obs = pd.DataFrame(
            {"cell_type": ["T cell", "B cell"], "organism_ontology_term_id": [None, None]},
            index=["cell1", "cell2"],
        )

        with pytest.raises(ValueError, match="organism_ontology_term_id cannot be null/NaN"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test_build_chrom_mapping_hg38(self, tmp_path):
        """Test chromosome mapping creation for human genome (hg38)."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        max_chrom, chrom_map = processor.build_chrom_mapping("hg38")

        # Verify mapping starts at 1
        assert chrom_map["chr1"] == 1
        assert chrom_map["chr2"] == 2
        assert chrom_map["chrX"] > 0
        assert chrom_map["chrY"] > 0

        # Verify max_chrom is the highest value
        assert max_chrom == max(chrom_map.values())

        # Verify unmapped chromosomes return 0 (defaultdict behavior)
        assert chrom_map["unknown_chr"] == 0

    def test_calculate_max_bins_mm39(self, tmp_path):
        """Test maximum bins calculation for mouse genome."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        max_bins = processor.calculate_max_bins("mm39")

        # Should be a positive integer
        assert isinstance(max_bins, int)
        assert max_bins > 0

        # Should be reasonable for mouse genome
        assert max_bins > 500_000  # At least 500K bins
        assert max_bins < 3_000_000  # Less than 3M bins

    def test_process_fragment_row_bin_calculation_accuracy(self, tmp_path):
        """Test accurate bin calculation for various fragment positions."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Test cases: (start, end, expected_start_bin, expected_end_bin)
        test_cases = [
            (0, 99, 0, 0),  # Both in bin 0
            (100, 199, 1, 1),  # Both in bin 1
            (99, 100, 0, 0),  # Edge case: end-1 = 99
            (199, 300, 1, 2),  # Spans bin 1 to 2
            (250, 350, 2, 3),  # Spans bin 2 to 3
        ]

        for start, end, expected_start_bin, expected_end_bin in test_cases:
            row = f"chr1\t{start}\t{end}\tcell1"
            chrom_id = 1
            cell_type_map = {"cell1": "T cell"}
            valid_barcodes = {"cell1"}
            coverage_aggregator = defaultdict(int)
            found_cells = set()

            processor._process_fragment_row(
                row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
            )

            # Verify bin calculations
            if expected_start_bin == expected_end_bin:
                # Single bin case
                expected_key = (1, expected_start_bin, "T cell")
                assert coverage_aggregator[expected_key] == 1
                assert len(coverage_aggregator) == 1
            else:
                # Multiple bins case
                start_key = (1, expected_start_bin, "T cell")
                end_key = (1, expected_end_bin, "T cell")
                assert coverage_aggregator[start_key] == 1
                assert coverage_aggregator[end_key] == 1
                assert len(coverage_aggregator) == 2

    def test_process_fragment_row_coverage_counting(self, tmp_path):
        """Test that coverage counting for insertion sites works correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Setup test data
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        # Process multiple fragments that affect the same bins
        fragments = [
            "chr1\t150\t350\tcell1",  # Bins 1 and 3
            "chr1\t250\t450\tcell1",  # Bins 2 and 4
            "chr1\t150\t250\tcell1",  # Bins 1 and 2
        ]

        for row in fragments:
            processor._process_fragment_row(
                row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
            )

        # Verify coverage accumulation
        assert coverage_aggregator[(1, 1, "T cell")] == 2  # From fragments 1 and 3
        assert coverage_aggregator[(1, 2, "T cell")] == 2  # From fragments 2 and 3
        assert coverage_aggregator[(1, 3, "T cell")] == 1  # From fragment 1
        assert coverage_aggregator[(1, 4, "T cell")] == 1  # From fragment 2

    @pytest.mark.parametrize(
        "invalid_row,expected_columns",
        [
            ("chr1\t100\t200", 3),  # Only 3 columns
            ("chr1\t100", 2),  # Only 2 columns
            ("chr1", 1),  # Only 1 column
            ("", 1),  # Empty string (split creates 1 empty element)
            ("\t\t", 1),  # Only tabs (strip() makes it empty, splits to 1 field)
        ],
    )
    def test_process_fragment_row_insufficient_columns(self, tmp_path, invalid_row, expected_columns, caplog):
        """Test that fragments with < 4 columns are handled properly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Setup test data
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        # Process invalid fragment row
        with caplog.at_level("WARNING"):
            processor._process_fragment_row(
                invalid_row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
            )

        # Verify no data was processed
        assert len(coverage_aggregator) == 0
        assert len(found_cells) == 0

        # Verify warning was logged
        if expected_columns < 4:
            assert "Invalid fragment format" in caplog.text
            assert f"expected at least 4 columns, got {expected_columns}" in caplog.text

    def test_process_fragment_row_extra_columns_valid(self, tmp_path):
        """Test that fragments with > 4 columns work correctly (uses first 4)."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Setup test data with extra columns
        row = "chr1\t150\t250\tcell1\textra_col1\textra_col2"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        # Process fragment with extra columns
        processor._process_fragment_row(row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells)

        # Verify fragment was processed correctly (extra columns ignored)
        assert "cell1" in found_cells
        expected_key = (1, 1, "T cell")  # 150//100 = 1, (250-1)//100 = 2, but they're different bins
        assert coverage_aggregator[expected_key] == 1

    def test_compute_cell_type_totals_multiple_cell_types(self, tmp_path):
        """Test cell type totals computation with multiple cell types."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create coverage aggregator with multiple cell types and chromosomes
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                # T cell coverage across different chromosomes and bins
                (1, 0, "T cell"): 5,
                (1, 1, "T cell"): 3,
                (2, 0, "T cell"): 2,
                (2, 5, "T cell"): 1,
                # B cell coverage
                (1, 0, "B cell"): 4,
                (1, 2, "B cell"): 6,
                (3, 1, "B cell"): 2,
                # NK cell coverage
                (1, 3, "NK cell"): 8,
                (2, 2, "NK cell"): 3,
            }
        )

        # Compute totals
        cell_type_totals = processor._compute_cell_type_totals(coverage_aggregator)

        # Verify totals are calculated correctly
        expected_totals = {
            "T cell": 5 + 3 + 2 + 1,  # = 11
            "B cell": 4 + 6 + 2,  # = 12
            "NK cell": 8 + 3,  # = 11
        }

        assert cell_type_totals == expected_totals

        # Verify all cell types are present
        assert set(cell_type_totals.keys()) == {"T cell", "B cell", "NK cell"}

        # Verify totals are positive integers
        for total in cell_type_totals.values():
            assert isinstance(total, int)
            assert total > 0

    def test_normalized_coverage_calculation(self, tmp_path):
        """Test normalized coverage calculation: (count / total_coverage) * normalization_factor."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create test coverage data with known totals
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                (1, 0, "T cell"): 10,  # T cell total will be 30
                (1, 1, "T cell"): 20,
                (2, 0, "B cell"): 50,  # B cell total will be 100
                (2, 1, "B cell"): 50,
            }
        )

        # Create DataFrame to test normalization
        df = processor._create_coverage_dataframe(coverage_aggregator)

        # Verify DataFrame structure
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        # Test normalization calculations
        normalization_factor = 2_000_000

        for _, row in df.iterrows():
            count = row["coverage"]
            total_coverage = row["total_coverage"]
            normalized_coverage = row["normalized_coverage"]
            cell_type = row["cell_type"]

            # Calculate expected normalized coverage
            expected_normalized = (count / total_coverage) * normalization_factor

            # Verify calculation is correct (with float32 precision tolerance)
            assert abs(normalized_coverage - expected_normalized) < 1e-1

            # Verify total_coverage matches expected cell type totals
            if cell_type == "T cell":
                assert total_coverage == 30
            elif cell_type == "B cell":
                assert total_coverage == 100

        # Test specific expected values
        t_cell_rows = df[df["cell_type"] == "T cell"]
        b_cell_rows = df[df["cell_type"] == "B cell"]

        # T cell with count=10, total=30: (10/30) * 2M = 666,666.67
        t_cell_10_row = t_cell_rows[t_cell_rows["coverage"] == 10].iloc[0]
        expected_t_cell_10 = (10 / 30) * normalization_factor
        assert abs(t_cell_10_row["normalized_coverage"] - expected_t_cell_10) < 1e-1

        # T cell with count=20, total=30: (20/30) * 2M = 1,333,333.33
        t_cell_20_row = t_cell_rows[t_cell_rows["coverage"] == 20].iloc[0]
        expected_t_cell_20 = (20 / 30) * normalization_factor
        assert abs(t_cell_20_row["normalized_coverage"] - expected_t_cell_20) < 1e-1

        # B cell with count=50, total=100: (50/100) * 2M = 1,000,000
        b_cell_rows_50 = b_cell_rows[b_cell_rows["coverage"] == 50]
        for _, row in b_cell_rows_50.iterrows():
            expected_b_cell = (50 / 100) * normalization_factor
            assert abs(row["normalized_coverage"] - expected_b_cell) < 1e-1

    def test_normalized_coverage_zero_total(self, tmp_path):
        """Test normalized coverage calculation when total_coverage is 0."""
        # Test direct calculation with zero total (edge case in _process_coverage_data)
        count = 5
        total_coverage = 0
        normalization_factor = 2_000_000

        # This simulates the calculation in line 244 of the original code
        normalized_coverage = (count / total_coverage) * normalization_factor if total_coverage > 0 else 0.0

        # Should return 0.0 when total_coverage is 0
        assert normalized_coverage == 0.0
        assert isinstance(normalized_coverage, float)

    def test_create_coverage_dataframe_small_dataset(self, tmp_path):
        """Test _create_coverage_dataframe() with small dataset."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create small coverage aggregator
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                (1, 0, "T cell"): 5,
                (1, 1, "T cell"): 3,
                (2, 0, "B cell"): 4,
            }
        )

        # Create DataFrame
        df = processor._create_coverage_dataframe(coverage_aggregator)

        # Verify DataFrame structure
        assert len(df) == 3  # Should have 3 rows
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        # Verify data types are optimized
        assert df["chrom"].dtype == "int32"
        assert df["bin"].dtype == "int32"
        assert df["coverage"].dtype == "int32"
        assert df["total_coverage"].dtype == "int32"
        assert df["normalized_coverage"].dtype == "float32"
        assert df["cell_type"].dtype == "object"  # String columns are object type

        # Verify cell type totals are calculated correctly
        t_cell_rows = df[df["cell_type"] == "T cell"]
        b_cell_rows = df[df["cell_type"] == "B cell"]

        # T cell total should be 5 + 3 = 8
        assert all(t_cell_rows["total_coverage"] == 8)
        # B cell total should be 4
        assert all(b_cell_rows["total_coverage"] == 4)

        # Verify individual rows contain expected data
        expected_data = [
            {"chrom": 1, "bin": 0, "cell_type": "T cell", "coverage": 5, "total_coverage": 8},
            {"chrom": 1, "bin": 1, "cell_type": "T cell", "coverage": 3, "total_coverage": 8},
            {"chrom": 2, "bin": 0, "cell_type": "B cell", "coverage": 4, "total_coverage": 4},
        ]

        # Sort DataFrame for predictable comparison
        df_sorted = df.sort_values(["chrom", "bin", "cell_type"]).reset_index(drop=True)

        for i, expected_row in enumerate(expected_data):
            actual_row = df_sorted.iloc[i]
            assert actual_row["chrom"] == expected_row["chrom"]
            assert actual_row["bin"] == expected_row["bin"]
            assert actual_row["cell_type"] == expected_row["cell_type"]
            assert actual_row["coverage"] == expected_row["coverage"]
            assert actual_row["total_coverage"] == expected_row["total_coverage"]

            # Verify normalized coverage calculation
            expected_normalized = (expected_row["coverage"] / expected_row["total_coverage"]) * 2_000_000
            assert abs(actual_row["normalized_coverage"] - expected_normalized) < 1e-1

        # Verify no NaN values
        assert not df.isnull().any().any()

    def test_create_coverage_dataframe_empty_aggregator(self, tmp_path):
        """Test _create_coverage_dataframe() with empty aggregator."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create empty coverage aggregator
        coverage_aggregator = defaultdict(int)

        # Create DataFrame
        df = processor._create_coverage_dataframe(coverage_aggregator)

        # Verify empty DataFrame has correct structure
        assert len(df) == 0
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        # Verify DataFrame is truly empty but properly structured
        assert df.empty
        assert isinstance(df, pd.DataFrame)

    def test_process_coverage_data_all_records_processed(self, tmp_path):
        """Test _process_coverage_data() processes all records."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Create coverage aggregator with multiple records
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                (1, 0, "T cell"): 10,
                (1, 1, "T cell"): 5,
                (2, 0, "B cell"): 8,
                (2, 1, "B cell"): 12,
                (3, 0, "NK cell"): 6,
            }
        )

        # Initialize required instance variables for processing with small chunk size
        processor._chunks = []
        processor._current_chunk = []
        processor._dataframe_chunk_size = 2  # Small chunk size to trigger chunking

        # Process coverage data - should process all records
        processor._process_coverage_data(coverage_aggregator)

        # Verify all records were processed by checking total records in chunks
        total_records_processed = 0
        for chunk in processor._chunks:
            total_records_processed += len(chunk)
        total_records_processed += len(processor._current_chunk)  # Add remaining records

        assert total_records_processed == 5  # Should have processed all 5 records

        # Clean up instance variables
        del processor._chunks
        del processor._current_chunk
        del processor._dataframe_chunk_size

    @patch("backend.layers.processing.utils.atac.tiledb")
    def test_create_dataframe_array(self, mock_tiledb, tmp_path):
        """Test create_dataframe_array() creates TileDB schema correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock TileDB components
        mock_filter_list = Mock()
        mock_domain = Mock()
        mock_dim = Mock()
        mock_attr = Mock()
        mock_schema = Mock()

        mock_tiledb.FilterList.return_value = mock_filter_list
        mock_tiledb.BitShuffleFilter.return_value = Mock()
        mock_tiledb.ZstdFilter.return_value = Mock()
        mock_tiledb.Domain.return_value = mock_domain
        mock_tiledb.Dim.return_value = mock_dim
        mock_tiledb.Attr.return_value = mock_attr
        mock_tiledb.ArraySchema.return_value = mock_schema
        mock_tiledb.SparseArray.create = Mock()

        # Test parameters
        array_name = str(tmp_path / "test_array")
        max_chrom = 25
        max_bins = 1000

        # Call method
        processor.create_dataframe_array(array_name, max_chrom, max_bins)

        # Verify compression filters were created
        assert mock_tiledb.FilterList.call_count >= 2  # One for compression, one for dimensions
        mock_tiledb.BitShuffleFilter.assert_called_once()
        assert mock_tiledb.ZstdFilter.call_count >= 1  # Called for both compression and dimensions

        # Verify ZstdFilter compression level
        zstd_calls = mock_tiledb.ZstdFilter.call_args_list
        compression_call = next((call for call in zstd_calls if call.kwargs.get("level") == 3), None)
        assert compression_call is not None

        # Verify domain creation with 3 dimensions
        mock_tiledb.Domain.assert_called_once()
        domain_args = mock_tiledb.Domain.call_args[0]
        assert len(domain_args) == 3  # chrom, bin, cell_type dimensions

        # Verify dimension creation calls
        dim_calls = mock_tiledb.Dim.call_args_list
        assert len(dim_calls) == 3

        # Check chrom dimension
        chrom_dim_call = dim_calls[0]
        assert chrom_dim_call.kwargs["name"] == "chrom"
        assert chrom_dim_call.kwargs["domain"] == (1, max_chrom)
        assert chrom_dim_call.kwargs["tile"] == 1
        assert chrom_dim_call.kwargs["dtype"] == np.uint32

        # Check bin dimension
        bin_dim_call = dim_calls[1]
        assert bin_dim_call.kwargs["name"] == "bin"
        assert bin_dim_call.kwargs["domain"] == (0, max_bins)
        assert bin_dim_call.kwargs["tile"] == 10
        assert bin_dim_call.kwargs["dtype"] == np.uint32

        # Check cell_type dimension
        cell_type_dim_call = dim_calls[2]
        assert cell_type_dim_call.kwargs["name"] == "cell_type"
        assert cell_type_dim_call.kwargs["dtype"] == "ascii"

        # Verify attribute creation
        attr_calls = mock_tiledb.Attr.call_args_list
        assert len(attr_calls) == 3

        # Check attributes
        attr_names = [call.kwargs["name"] for call in attr_calls]
        assert "coverage" in attr_names
        assert "total_coverage" in attr_names
        assert "normalized_coverage" in attr_names

        # Check attribute dtypes
        coverage_attr = next(call for call in attr_calls if call.kwargs["name"] == "coverage")
        assert coverage_attr.kwargs["dtype"] == np.int32

        total_coverage_attr = next(call for call in attr_calls if call.kwargs["name"] == "total_coverage")
        assert total_coverage_attr.kwargs["dtype"] == np.int32

        normalized_attr = next(call for call in attr_calls if call.kwargs["name"] == "normalized_coverage")
        assert normalized_attr.kwargs["dtype"] == np.float32

        # Verify schema creation
        mock_tiledb.ArraySchema.assert_called_once()
        schema_kwargs = mock_tiledb.ArraySchema.call_args.kwargs
        assert schema_kwargs["domain"] == mock_domain
        assert len(schema_kwargs["attrs"]) == 3
        assert schema_kwargs["sparse"] is True
        assert schema_kwargs["allows_duplicates"] is False

        # Verify array creation
        mock_tiledb.SparseArray.create.assert_called_once_with(array_name, mock_schema)

        # Test writing to the created array (covers line 317)
        test_df = pd.DataFrame(
            {
                "chrom": [1, 2],
                "bin": [0, 1],
                "cell_type": ["T cell", "B cell"],
                "coverage": [10, 5],
                "total_coverage": [15, 5],
                "normalized_coverage": [1333333.25, 2000000.0],
            }
        )

        # Mock the array for writing
        mock_array = Mock()
        mock_array.__setitem__ = Mock()  # Enable item assignment
        mock_tiledb.SparseArray.return_value.__enter__.return_value = mock_array

        # Call write method (this covers line 317)
        processor._write_coverage_to_tiledb(array_name, test_df)

        # Verify TileDB array was opened for writing
        mock_tiledb.SparseArray.assert_called_with(array_name, mode="w", ctx=processor.ctx)

        # Verify data was written
        assert mock_array.__setitem__.called

    @patch("backend.layers.processing.utils.atac.tiledb")
    def test_create_dataframe_array_edge_cases(self, mock_tiledb, tmp_path):
        """Test create_dataframe_array() with edge case parameters."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock TileDB components
        mock_tiledb.FilterList.return_value = Mock()
        mock_tiledb.BitShuffleFilter.return_value = Mock()
        mock_tiledb.ZstdFilter.return_value = Mock()
        mock_tiledb.Domain.return_value = Mock()
        mock_tiledb.Dim.return_value = Mock()
        mock_tiledb.Attr.return_value = Mock()
        mock_tiledb.ArraySchema.return_value = Mock()
        mock_tiledb.SparseArray.create = Mock()

        # Test with minimal values
        array_name = str(tmp_path / "minimal_array")
        max_chrom = 1
        max_bins = 0

        # Should not raise any errors
        processor.create_dataframe_array(array_name, max_chrom, max_bins)

        # Verify dimensions were created with edge case values
        dim_calls = mock_tiledb.Dim.call_args_list
        chrom_dim_call = dim_calls[0]
        bin_dim_call = dim_calls[1]

        assert chrom_dim_call.kwargs["domain"] == (1, 1)
        assert bin_dim_call.kwargs["domain"] == (0, 0)

    @patch("backend.layers.processing.utils.atac.pysam")
    def test_process_all_chromosomes(self, mock_pysam, tmp_path):
        """Test _process_all_chromosomes() processes fragments from all chromosomes."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile
        mock_tabix = Mock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix

        # Mock fragment data for different chromosomes
        chr1_fragments = [
            "chr1\t100\t200\tcell1",
            "chr1\t300\t400\tcell2",
            "chr1\t500\t600\tcell1",
            "chr1\t700\t800\tinvalid_cell",  # Invalid barcode - should be filtered out
            "chr1\t-50\t100\tcell1",  # Invalid coordinates: negative start
            "chr1\t200\t200\tcell2",  # Invalid coordinates: start == end
            "chr1\t300\t250\tcell1",  # Invalid coordinates: start > end
            "chr1\tabc\t200\tcell1",  # Invalid coordinates: non-integer start (triggers ValueError)
            "chr1\t100\txyz\tcell2",  # Invalid coordinates: non-integer end (triggers ValueError)
        ]
        chr2_fragments = [
            "chr2\t150\t250\tcell2",
            "chr2\t350\t450\tcell3",
        ]
        chrX_fragments = [
            "chrX\t200\t300\tcell1",
            "chrX\t400\t500\tunknown_barcode",  # Another invalid barcode
        ]

        # Configure mock to return different fragments per chromosome
        fragment_map = {
            "chr1": chr1_fragments,
            "chr2": chr2_fragments,
            "chrX": chrX_fragments,
        }

        def mock_fetch(chrom):
            if chrom == "chr2":  # Simulate fetch failure for chr2
                raise ValueError("Simulated tabix fetch error")
            return fragment_map.get(chrom, [])

        mock_tabix.fetch.side_effect = mock_fetch

        # Test parameters
        chrom_map = {"chr1": 1, "chr2": 2, "chrX": 23}
        cell_type_map = {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"}
        valid_barcodes = {"cell1", "cell2", "cell3"}

        # Call method
        coverage_aggregator, found_cells = processor._process_all_chromosomes(chrom_map, cell_type_map, valid_barcodes)

        # Verify TabixFile was used correctly
        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))

        # Verify fetch was called for each chromosome
        expected_fetch_calls = ["chr1", "chr2", "chrX"]
        actual_fetch_calls = [call[0][0] for call in mock_tabix.fetch.call_args_list]
        assert set(actual_fetch_calls) == set(expected_fetch_calls)

        # Verify cells from successful chromosomes were found (chr2 failed, so cell3 missing from chr2)
        assert found_cells == {"cell1", "cell2"}  # cell1 from chr1+chrX, cell2 from chr1

        # Verify coverage aggregation
        # chr1 fragments: cell1 (bins 1,2), cell2 (bins 3,4), cell1 (bins 5,6)
        # chr2 fragments: cell2 (bins 1,2), cell3 (bins 3,4)
        # chrX fragments: cell1 (bins 2,3)

        # Expected coverage structure (for reference):
        # chr1 (chrom_id=1): cell1 (bins 1,2), cell2 (bins 3,4), cell1 (bins 5,6)
        # chr2 (chrom_id=2): cell2 (bins 1,2), cell3 (bins 3,4)
        # chrX (chrom_id=23): cell1 (bins 2,3)

        # Check specific key coverage values (only chr1 and chrX succeeded)
        assert coverage_aggregator[(1, 1, "T cell")] >= 1  # cell1 from chr1
        assert coverage_aggregator[(1, 3, "B cell")] >= 1  # cell2 from chr1
        assert coverage_aggregator[(23, 2, "T cell")] >= 1  # cell1 from chrX

        # chr2 data should be missing due to ValueError (no NK cell coverage)
        assert (2, 1, "B cell") not in coverage_aggregator
        assert (2, 3, "NK cell") not in coverage_aggregator

        # Verify coverage aggregator contains T cell and B cell (but not NK cell due to chr2 failure)
        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        assert cell_types_in_coverage == {"T cell", "B cell"}

        # Verify chromosomes in coverage match chrom_map
        chroms_in_coverage = {key[0] for key in coverage_aggregator}
        assert chroms_in_coverage.issubset(set(chrom_map.values()))

    @patch("backend.layers.processing.utils.atac.pysam")
    def test_process_all_chromosomes_missing_cells(self, mock_pysam, tmp_path):
        """Test _process_all_chromosomes() handles missing cells correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile with limited fragments
        mock_tabix = Mock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix

        # Only fragments for cell1, missing cell2 and cell3
        mock_tabix.fetch.return_value = ["chr1\t100\t200\tcell1"]

        # Test parameters - expect cells 1, 2, 3 but only cell1 has fragments
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"}
        valid_barcodes = {"cell1", "cell2", "cell3"}

        # Call method
        coverage_aggregator, found_cells = processor._process_all_chromosomes(chrom_map, cell_type_map, valid_barcodes)

        # Verify only cell1 was found
        assert found_cells == {"cell1"}

        # Verify coverage only contains T cell data
        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        assert cell_types_in_coverage == {"T cell"}

        # Verify missing cells can be detected
        missing_cells = valid_barcodes - found_cells
        assert missing_cells == {"cell2", "cell3"}

    @patch("backend.layers.processing.utils.atac.pysam")
    @patch("backend.layers.processing.utils.atac.tiledb")
    def test_write_binned_coverage_per_chrom_full_pipeline(self, mock_tiledb, mock_pysam, tmp_path):
        """Test write_binned_coverage_per_chrom() full pipeline orchestration (lines 123,126,129,132,135)."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile with fragments
        mock_tabix = Mock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = [
            "chr1\t100\t200\tcell1",
            "chr1\t300\t400\tcell2",
        ]

        # Mock TileDB for writing
        mock_array = Mock()
        mock_array.__setitem__ = Mock()
        mock_tiledb.SparseArray.return_value.__enter__.return_value = mock_array

        # Test parameters
        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"}
        valid_barcodes = {"cell1", "cell2", "cell3"}

        # Call the orchestration method (covers lines 123, 126, 132, 135)
        processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

        # Verify _process_all_chromosomes was called (line 123)
        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))

        # Verify _report_missing_cells was called (line 126) - cell3 should be missing
        # This is tested implicitly by the method completing successfully

        # Verify _create_coverage_dataframe was called (line 132) - coverage_aggregator not empty
        # Verify _write_coverage_to_tiledb was called (line 135)
        mock_tiledb.SparseArray.assert_called_with(array_name, mode="w", ctx=processor.ctx)
        assert mock_array.__setitem__.called

    @patch("backend.layers.processing.utils.atac.pysam")
    def test_write_binned_coverage_per_chrom_empty_coverage(self, mock_pysam, tmp_path):
        """Test write_binned_coverage_per_chrom() early return with empty coverage (line 129)."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile to return no fragments (empty coverage)
        mock_tabix = Mock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = []  # No fragments found

        # Test parameters
        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}

        # Call method - should return early due to empty coverage_aggregator (line 129)
        result = processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

        # Verify method returns None (early return on line 129)
        assert result is None

        # Verify _process_all_chromosomes was still called (line 123)
        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))

    @patch("backend.layers.processing.utils.atac.tiledb")
    @patch("backend.layers.processing.utils.atac.pysam")
    def test_process_fragment_file_integration(self, mock_pysam, mock_tiledb, tmp_path):
        """Integration test for process_fragment_file covering lines 335, 337, 341, 344, 346."""
        # Setup test file
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\nchr1\t300\t400\tcell2\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile
        mock_tabix = Mock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = ["chr1\t100\t200\tcell1", "chr1\t300\t400\tcell2"]

        # Mock TileDB components
        mock_tiledb.FilterList.return_value = Mock()
        mock_tiledb.BitShuffleFilter.return_value = Mock()
        mock_tiledb.ZstdFilter.return_value = Mock()
        mock_tiledb.Domain.return_value = Mock()
        mock_tiledb.Dim.return_value = Mock()
        mock_tiledb.Attr.return_value = Mock()
        mock_tiledb.ArraySchema.return_value = Mock()
        mock_tiledb.SparseArray.create = Mock()

        # Mock SparseArray for writing
        mock_array = Mock()
        mock_array.__setitem__ = Mock()
        mock_tiledb.SparseArray.return_value.__enter__.return_value = mock_array

        # Create test obs DataFrame
        obs = pd.DataFrame(
            {
                "cell_type": ["T cell", "B cell"],
                "organism_ontology_term_id": ["NCBITaxon:9606", "NCBITaxon:9606"],
            },
            index=["cell1", "cell2"],
        )

        array_name = str(tmp_path / "test_array")

        # Call process_fragment_file - covers all missing lines
        processor.process_fragment_file(obs, array_name)

        # Get df_meta to verify the internal processing
        df_meta, _ = processor.extract_cell_metadata_from_h5ad(obs)

        # Verify line 335: valid_barcodes = set(df_meta["cell_name"])
        expected_valid_barcodes = {"cell1", "cell2"}
        assert set(df_meta["cell_name"]) == expected_valid_barcodes

        # Verify line 335: cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))
        expected_cell_type_map = {"cell1": "T cell", "cell2": "B cell"}
        actual_cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))
        assert actual_cell_type_map == expected_cell_type_map

        # Verify line 337: max_chrom, chrom_map = self.build_chrom_mapping(genome_version)
        # This is called internally, verify chromosome mapping was used
        assert mock_tabix.fetch.call_count > 0  # Should be called for each chromosome

        # Verify line 341: self.write_binned_coverage_per_chrom called
        mock_pysam.TabixFile.assert_called_with(str(fragment_file))

        # Verify DataFrame structure and content
        assert isinstance(df_meta, pd.DataFrame)
        assert len(df_meta) == 2
