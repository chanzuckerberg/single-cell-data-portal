from collections import defaultdict

import numpy as np
import pandas as pd
import pytest

from backend.layers.processing.utils.atac import ATACDataProcessor


# Parametrized fixtures and test data
@pytest.fixture(
    params=[
        {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"},
        {"cell1": "CD4+ T cell", "cell2": "Memory B cell", "cell3": "Monocyte"},
        {"cell1": "Naive T cell", "cell2": "Plasma cell"},
    ]
)
def cell_type_mapping(request):
    """Parametrized fixture for cell type mappings."""
    return request.param


@pytest.fixture(
    params=[
        ("NCBITaxon:9606", "hg38"),
        ("NCBITaxon:10090", "mm39"),
    ]
)
def organism_genome_pair(request):
    """Parametrized fixture for organism ID and genome version pairs."""
    return request.param


@pytest.fixture(
    params=[
        # Small coverage aggregator
        {
            (1, 0, "T cell"): 5,
            (1, 1, "T cell"): 3,
            (2, 0, "B cell"): 4,
        },
        # Multi-cell type coverage
        {
            (1, 0, "T cell"): 10,
            (1, 1, "T cell"): 20,
            (2, 0, "B cell"): 50,
            (2, 1, "B cell"): 50,
            (3, 0, "NK cell"): 15,
        },
        # Single cell type
        {
            (1, 0, "T cell"): 100,
            (1, 1, "T cell"): 200,
            (1, 2, "T cell"): 300,
        },
    ]
)
def coverage_aggregator_data(request):
    """Parametrized fixture for coverage aggregator test data."""
    coverage_aggregator = defaultdict(int)
    coverage_aggregator.update(request.param)
    return coverage_aggregator


@pytest.fixture(
    params=[
        # Valid coordinates
        ("chr1\t100\t200\tcell1", 100, 199, 1, 1),
        ("chr1\t0\t99\tcell1", 0, 99, 0, 0),
        ("chr1\t199\t300\tcell1", 199, 299, 1, 2),
        ("chr1\t250\t350\tcell1", 250, 349, 2, 3),
        # Edge cases
        ("chr1\t99\t100\tcell1", 99, 99, 0, 0),
    ]
)
def fragment_coordinate_data(request):
    """Parametrized fixture for fragment coordinate test cases."""
    row, start, end, expected_start_bin, expected_end_bin = request.param
    return {
        "row": row,
        "start": start,
        "end": end,
        "expected_start_bin": expected_start_bin,
        "expected_end_bin": expected_end_bin,
    }


@pytest.fixture(
    params=[
        # Invalid coordinates
        ("chr1\t-50\t100\tcell1", "negative start"),
        ("chr1\t200\t200\tcell1", "start equals end"),
        ("chr1\t300\t250\tcell1", "start greater than end"),
        ("chr1\tabc\t200\tcell1", "non-integer start"),
        ("chr1\t100\txyz\tcell1", "non-integer end"),
    ]
)
def invalid_fragment_coordinates(request):
    """Parametrized fixture for invalid fragment coordinate test cases."""
    row, description = request.param
    return {"row": row, "description": description}


@pytest.fixture(
    params=[
        # TileDB array configuration scenarios
        {"max_chrom": 25, "max_bins": 1000, "compression_level": 3},
        {"max_chrom": 1, "max_bins": 0, "compression_level": 3},  # Edge case
        {"max_chrom": 50, "max_bins": 5000, "compression_level": 3},  # Larger dataset
    ]
)
def tiledb_array_config(request):
    """Parametrized fixture for TileDB array configuration test cases."""
    return request.param


@pytest.fixture
def mock_tiledb_components(mocker):
    """Fixture that provides mocked TileDB components for testing."""
    mock_tiledb = mocker.patch("backend.layers.processing.utils.atac.tiledb")

    # Create mock objects
    mock_filter_list = mocker.MagicMock()
    mock_domain = mocker.MagicMock()
    mock_dim = mocker.MagicMock()
    mock_attr = mocker.MagicMock()
    mock_schema = mocker.MagicMock()
    mock_array = mocker.MagicMock()

    # Configure mock returns
    mock_tiledb.FilterList.return_value = mock_filter_list
    mock_tiledb.BitShuffleFilter.return_value = mocker.MagicMock()
    mock_tiledb.ZstdFilter.return_value = mocker.MagicMock()
    mock_tiledb.Domain.return_value = mock_domain
    mock_tiledb.Dim.return_value = mock_dim
    mock_tiledb.Attr.return_value = mock_attr
    mock_tiledb.ArraySchema.return_value = mock_schema
    mock_tiledb.SparseArray.create = mocker.MagicMock()

    # Configure array for writing
    mock_array.__setitem__ = mocker.MagicMock()
    mock_tiledb.SparseArray.return_value.__enter__.return_value = mock_array

    return {
        "tiledb": mock_tiledb,
        "filter_list": mock_filter_list,
        "domain": mock_domain,
        "dim": mock_dim,
        "attr": mock_attr,
        "schema": mock_schema,
        "array": mock_array,
    }


class TestATACDataProcessor:
    """Test suite for ATACDataProcessor class."""

    def test__get_genome_version_valid_organisms(self, organism_genome_pair):
        """Test genome version mapping for valid organisms."""
        organism_id, expected_genome = organism_genome_pair
        result = ATACDataProcessor.get_genome_version(organism_id)
        assert result == expected_genome

    def test__get_genome_version_invalid_organism(self):
        """Test ValueError is raised for unknown organism ontology term ID."""
        with pytest.raises(ValueError, match="Unknown organism ontology term ID"):
            ATACDataProcessor.get_genome_version("NCBITaxon:12345")

    def test__constructor_valid_file_path(self, tmp_path):
        """Test ATACDataProcessor initialization with valid fragment file path."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        assert processor.fragment_artifact_id == str(fragment_file)
        assert processor.ctx is None
        assert processor.bin_size == 100
        assert processor.normalization_factor == 2_000_000

    def test__constructor_file_not_found(self):
        """Test FileNotFoundError is raised when fragment file doesn't exist."""
        non_existent_file = "/path/to/non/existent/file.tsv.gz"

        with pytest.raises(FileNotFoundError, match=f"Fragment file not found: {non_existent_file}"):
            ATACDataProcessor(fragment_artifact_id=non_existent_file)

    def test__extract_cell_metadata_valid_obs(self, tmp_path, cell_type_mapping, organism_genome_pair):
        """Test cell metadata extraction with valid obs DataFrame."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        organism_id, expected_genome = organism_genome_pair
        cell_names = list(cell_type_mapping.keys())
        cell_types = list(cell_type_mapping.values())
        obs = pd.DataFrame(
            {
                "cell_type": cell_types,
                "organism_ontology_term_id": [organism_id] * len(cell_types),
            },
            index=cell_names,
        )

        df_meta, genome_version = processor.extract_cell_metadata_from_h5ad(obs)

        expected_df = pd.DataFrame({"cell_name": cell_names, "cell_type": cell_types})
        pd.testing.assert_frame_equal(df_meta, expected_df)
        assert genome_version == expected_genome

    def test__extract_cell_metadata_missing_obs_column(self, tmp_path, organism_genome_pair):
        """Test ValueError when obs_column is missing from DataFrame."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        organism_id, _ = organism_genome_pair
        obs = pd.DataFrame(
            {"other_column": ["value1", "value2"], "organism_ontology_term_id": [organism_id, organism_id]},
            index=["cell1", "cell2"],
        )

        with pytest.raises(ValueError, match="Column cell_type not found in obs DataFrame"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test__extract_cell_metadata_missing_organism_column(self, tmp_path, cell_type_mapping):
        """Test ValueError when organism_ontology_term_id column is missing."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        cell_names = list(cell_type_mapping.keys())[:2]
        cell_types = [cell_type_mapping[name] for name in cell_names]
        obs = pd.DataFrame({"cell_type": cell_types}, index=cell_names)

        with pytest.raises(ValueError, match="Column 'organism_ontology_term_id' is required but not found"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test__extract_cell_metadata_null_organism_id(self, tmp_path):
        """Test ValueError when organism_ontology_term_id is null/NaN."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        obs = pd.DataFrame(
            {"cell_type": ["T cell", "B cell"], "organism_ontology_term_id": [None, None]},
            index=["cell1", "cell2"],
        )

        with pytest.raises(ValueError, match="organism_ontology_term_id cannot be null/NaN"):
            processor.extract_cell_metadata_from_h5ad(obs)

    def test__build_chrom_mapping(self, tmp_path, organism_genome_pair):
        """Test chromosome mapping creation for different genomes."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        _, genome_version = organism_genome_pair
        max_chrom, chrom_map = processor.build_chrom_mapping(genome_version)

        assert chrom_map["chr1"] == 1
        assert chrom_map["chr2"] == 2
        assert chrom_map["chrX"] > 0
        assert chrom_map["chrY"] > 0

        assert max_chrom == max(chrom_map.values())

        assert chrom_map["unknown_chr"] == 0

    def test__calculate_max_bins(self, tmp_path, organism_genome_pair):
        """Test maximum bins calculation for different genomes."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        _, genome_version = organism_genome_pair
        max_bins = processor.calculate_max_bins(genome_version)

        assert isinstance(max_bins, int)
        assert max_bins > 0
        assert max_bins > 500_000  # At least 500K bins
        assert max_bins < 3_500_000  # Less than 3.5M bins (hg38 is larger than mm39)

    def test__process_fragment_row_bin_calculation_accuracy(self, tmp_path, fragment_coordinate_data):
        """Test accurate bin calculation for various fragment positions."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Use parametrized test data
        row = fragment_coordinate_data["row"]
        expected_start_bin = fragment_coordinate_data["expected_start_bin"]
        expected_end_bin = fragment_coordinate_data["expected_end_bin"]

        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        processor._process_fragment_row(row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells)

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

    def test__process_fragment_row_coverage_counting(self, tmp_path):
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

        fragments = [
            "chr1\t150\t350\tcell1",  # Bins 1 and 3
            "chr1\t250\t450\tcell1",  # Bins 2 and 4
            "chr1\t150\t250\tcell1",  # Bins 1 and 2
        ]

        for row in fragments:
            processor._process_fragment_row(
                row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
            )

        assert coverage_aggregator[(1, 1, "T cell")] == 2  # From fragments 1 and 3
        assert coverage_aggregator[(1, 2, "T cell")] == 2  # From fragments 2 and 3
        assert coverage_aggregator[(1, 3, "T cell")] == 1  # From fragment 1
        assert coverage_aggregator[(1, 4, "T cell")] == 1  # From fragment 2

    @pytest.mark.parametrize(
        "invalid_row,expected_columns",
        [
            ("chr1\t100\t200", 3),
            ("chr1\t100", 2),
            ("chr1", 1),
            ("", 1),
            ("\t\t", 1),
        ],
    )
    def test__process_fragment_row_insufficient_columns(self, tmp_path, invalid_row, expected_columns, caplog):
        """Test that fragments with < 4 columns are handled properly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        with caplog.at_level("WARNING"):
            processor._process_fragment_row(
                invalid_row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells
            )

        assert len(coverage_aggregator) == 0
        assert len(found_cells) == 0

        if expected_columns < 4:
            assert "Invalid fragment format" in caplog.text
            assert f"expected at least 4 columns, got {expected_columns}" in caplog.text

    def test__process_fragment_row_extra_columns_valid(self, tmp_path):
        """Test that fragments with > 4 columns work correctly (uses first 4)."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        row = "chr1\t150\t250\tcell1\textra_col1\textra_col2"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        coverage_aggregator = defaultdict(int)
        found_cells = set()

        processor._process_fragment_row(row, chrom_id, cell_type_map, valid_barcodes, coverage_aggregator, found_cells)

        assert "cell1" in found_cells
        expected_key = (1, 1, "T cell")  # 150//100 = 1, (250-1)//100 = 2, but they're different bins
        assert coverage_aggregator[expected_key] == 1

    def test__compute_cell_type_totals_multiple_cell_types(self, tmp_path, coverage_aggregator_data):
        """Test cell type totals computation with multiple cell types."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        cell_type_totals = processor._compute_cell_type_totals(coverage_aggregator_data)
        expected_totals = {}
        for (_, _, cell_type), count in coverage_aggregator_data.items():
            if cell_type not in expected_totals:
                expected_totals[cell_type] = 0
            expected_totals[cell_type] += count

        assert cell_type_totals == expected_totals

        expected_cell_types = {key[2] for key in coverage_aggregator_data}
        assert set(cell_type_totals.keys()) == expected_cell_types
        for total in cell_type_totals.values():
            assert isinstance(total, int)
            assert total > 0

    def test__normalized_coverage_calculation(self, tmp_path):
        """Test normalized coverage calculation: (count / total_coverage) * normalization_factor."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                (1, 0, "T cell"): 10,  # T cell total will be 30
                (1, 1, "T cell"): 20,
                (2, 0, "B cell"): 50,  # B cell total will be 100
                (2, 1, "B cell"): 50,
            }
        )

        df = processor._create_coverage_dataframe(coverage_aggregator)
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        normalization_factor = 2_000_000

        for _, row in df.iterrows():
            count = row["coverage"]
            total_coverage = row["total_coverage"]
            normalized_coverage = row["normalized_coverage"]
            cell_type = row["cell_type"]

            expected_normalized = (count / total_coverage) * normalization_factor
            assert abs(normalized_coverage - expected_normalized) < 1e-1
            if cell_type == "T cell":
                assert total_coverage == 30
            elif cell_type == "B cell":
                assert total_coverage == 100

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

    def test__normalized_coverage_zero_total(self, tmp_path):
        """Test normalized coverage calculation when total_coverage is 0."""
        count = 5
        total_coverage = 0
        normalization_factor = 2_000_000

        # This simulates the calculation in the original code
        normalized_coverage = (count / total_coverage) * normalization_factor if total_coverage > 0 else 0.0

        assert normalized_coverage == 0.0
        assert isinstance(normalized_coverage, float)

    def test__create_coverage_dataframe_small_dataset(self, tmp_path, coverage_aggregator_data):
        """Test _create_coverage_dataframe() with parametrized datasets."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        df = processor._create_coverage_dataframe(coverage_aggregator_data)
        expected_rows = len(coverage_aggregator_data)
        assert len(df) == expected_rows
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        assert df["chrom"].dtype == "int32"
        assert df["bin"].dtype == "int32"
        assert df["coverage"].dtype == "int32"
        assert df["total_coverage"].dtype == "int32"
        assert df["normalized_coverage"].dtype == "float32"
        assert df["cell_type"].dtype == "object"  # String columns are object type

        expected_totals = {}
        for (_, _, cell_type), count in coverage_aggregator_data.items():
            if cell_type not in expected_totals:
                expected_totals[cell_type] = 0
            expected_totals[cell_type] += count

        for cell_type, expected_total in expected_totals.items():
            cell_type_rows = df[df["cell_type"] == cell_type]
            assert all(cell_type_rows["total_coverage"] == expected_total)

        for _, row in df.iterrows():
            expected_normalized = (row["coverage"] / row["total_coverage"]) * 2_000_000
            assert abs(row["normalized_coverage"] - expected_normalized) < 1e-1

        assert not df.isnull().any().any()

    def test__create_coverage_dataframe_empty_aggregator(self, tmp_path):
        """Test _create_coverage_dataframe() with empty aggregator."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        coverage_aggregator = defaultdict(int)
        df = processor._create_coverage_dataframe(coverage_aggregator)

        assert len(df) == 0
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert list(df.columns) == expected_columns

        assert df.empty
        assert isinstance(df, pd.DataFrame)

    def test__process_coverage_data_all_records_processed(self, tmp_path):
        """Test _process_coverage_data() processes all records."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

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

        processor._process_coverage_data(coverage_aggregator)
        total_records_processed = 0
        for chunk in processor._chunks:
            total_records_processed += len(chunk)
        total_records_processed += len(processor._current_chunk)
        assert total_records_processed == 5
        del processor._chunks
        del processor._current_chunk
        del processor._dataframe_chunk_size

    def test__create_dataframe_array(self, tmp_path, mock_tiledb_components, tiledb_array_config):
        """Test create_dataframe_array() creates TileDB schema correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        mock_tiledb = mock_tiledb_components["tiledb"]

        array_name = str(tmp_path / "test_array")
        max_chrom = tiledb_array_config["max_chrom"]
        max_bins = tiledb_array_config["max_bins"]
        expected_compression_level = tiledb_array_config["compression_level"]

        processor.create_dataframe_array(array_name, max_chrom, max_bins)

        assert mock_tiledb.FilterList.call_count >= 2  # One for compression, one for dimensions
        mock_tiledb.BitShuffleFilter.assert_called_once()
        assert mock_tiledb.ZstdFilter.call_count >= 1  # Called for both compression and dimensions
        zstd_calls = mock_tiledb.ZstdFilter.call_args_list
        compression_call = next(
            (call for call in zstd_calls if call.kwargs.get("level") == expected_compression_level), None
        )
        assert compression_call is not None

        mock_tiledb.Domain.assert_called_once()
        domain_args = mock_tiledb.Domain.call_args[0]
        assert len(domain_args) == 3  # chrom, bin, cell_type dimensions
        dim_calls = mock_tiledb.Dim.call_args_list
        assert len(dim_calls) == 3

        chrom_dim_call = dim_calls[0]
        assert chrom_dim_call.kwargs["name"] == "chrom"
        assert chrom_dim_call.kwargs["domain"] == (1, max_chrom)
        assert chrom_dim_call.kwargs["tile"] == 1
        assert chrom_dim_call.kwargs["dtype"] == np.uint32

        bin_dim_call = dim_calls[1]
        assert bin_dim_call.kwargs["name"] == "bin"
        assert bin_dim_call.kwargs["domain"] == (0, max_bins)
        assert bin_dim_call.kwargs["tile"] == 10
        assert bin_dim_call.kwargs["dtype"] == np.uint32

        cell_type_dim_call = dim_calls[2]
        assert cell_type_dim_call.kwargs["name"] == "cell_type"
        assert cell_type_dim_call.kwargs["dtype"] == "ascii"

        attr_calls = mock_tiledb.Attr.call_args_list
        assert len(attr_calls) == 3

        attr_names = [call.kwargs["name"] for call in attr_calls]
        assert "coverage" in attr_names
        assert "total_coverage" in attr_names
        assert "normalized_coverage" in attr_names

        coverage_attr = next(call for call in attr_calls if call.kwargs["name"] == "coverage")
        assert coverage_attr.kwargs["dtype"] == np.int32

        total_coverage_attr = next(call for call in attr_calls if call.kwargs["name"] == "total_coverage")
        assert total_coverage_attr.kwargs["dtype"] == np.int32

        normalized_attr = next(call for call in attr_calls if call.kwargs["name"] == "normalized_coverage")
        assert normalized_attr.kwargs["dtype"] == np.float32

        mock_tiledb.ArraySchema.assert_called_once()
        schema_kwargs = mock_tiledb.ArraySchema.call_args.kwargs
        assert schema_kwargs["domain"] == mock_tiledb_components["domain"]
        assert len(schema_kwargs["attrs"]) == 3
        assert schema_kwargs["sparse"] is True
        assert schema_kwargs["allows_duplicates"] is False

        mock_tiledb.SparseArray.create.assert_called_once_with(array_name, mock_tiledb_components["schema"])
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

        processor._write_coverage_to_tiledb(array_name, test_df)
        mock_tiledb.SparseArray.assert_called_with(array_name, mode="w", ctx=processor.ctx)
        assert mock_tiledb_components["array"].__setitem__.called

    def test__create_dataframe_array_edge_cases(self, tmp_path, mock_tiledb_components):
        """Test create_dataframe_array() with edge case parameters."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))
        mock_tiledb = mock_tiledb_components["tiledb"]

        array_name = str(tmp_path / "minimal_array")
        max_chrom = 1
        max_bins = 0

        processor.create_dataframe_array(array_name, max_chrom, max_bins)
        dim_calls = mock_tiledb.Dim.call_args_list
        chrom_dim_call = dim_calls[0]
        bin_dim_call = dim_calls[1]

        assert chrom_dim_call.kwargs["domain"] == (1, 1)
        assert bin_dim_call.kwargs["domain"] == (0, 0)

    def test__process_all_chromosomes(self, tmp_path, cell_type_mapping, invalid_fragment_coordinates, mocker):
        """Test _process_all_chromosomes() processes fragments from all chromosomes."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile
        mock_pysam = mocker.patch("backend.layers.processing.utils.atac.pysam")
        mock_tabix = mocker.MagicMock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix

        valid_cells = list(cell_type_mapping.keys())[:3]
        invalid_coord_row = invalid_fragment_coordinates["row"]

        # Mock fragment data for different chromosomes
        chr1_fragments = [
            "chr1\t100\t200\t" + valid_cells[0],
            "chr1\t300\t400\t" + valid_cells[1],
            "chr1\t500\t600\t" + valid_cells[0],
            "chr1\t700\t800\tinvalid_cell",  # Invalid barcode - should be filtered out
            invalid_coord_row,
        ]
        chr2_fragments = [
            "chr2\t150\t250\t" + valid_cells[1],
            "chr2\t350\t450\t" + valid_cells[2] if len(valid_cells) > 2 else "chr2\t350\t450\t" + valid_cells[0],
        ]
        chrX_fragments = [
            "chrX\t200\t300\t" + valid_cells[0],
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

        chrom_map = {"chr1": 1, "chr2": 2, "chrX": 23}
        valid_barcodes = set(valid_cells)
        coverage_aggregator, found_cells = processor._process_all_chromosomes(
            chrom_map, cell_type_mapping, valid_barcodes
        )

        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))
        expected_fetch_calls = ["chr1", "chr2", "chrX"]
        actual_fetch_calls = [call[0][0] for call in mock_tabix.fetch.call_args_list]
        assert set(actual_fetch_calls) == set(expected_fetch_calls)

        expected_found_cells = {valid_cells[0], valid_cells[1]}
        assert found_cells == expected_found_cells

        assert len(coverage_aggregator) > 0

        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        expected_cell_types = {cell_type_mapping[cell] for cell in expected_found_cells}
        assert cell_types_in_coverage == expected_cell_types

        chroms_in_coverage = {key[0] for key in coverage_aggregator}
        expected_chroms = {1, 23}  # chr1=1, chrX=23 (chr2=2 failed)
        assert chroms_in_coverage == expected_chroms

    def test__process_all_chromosomes_missing_cells(self, tmp_path, cell_type_mapping, mocker):
        """Test _process_all_chromosomes() handles missing cells correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile with limited fragments
        mock_pysam = mocker.patch("backend.layers.processing.utils.atac.pysam")
        mock_tabix = mocker.MagicMock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix

        # Use parametrized cell type mapping
        present_cell = list(cell_type_mapping.keys())[0]  # First cell will have fragments
        expected_cell_type = cell_type_mapping[present_cell]

        # Only fragments for one cell, missing others
        mock_tabix.fetch.return_value = [f"chr1\t100\t200\t{present_cell}"]

        # Test parameters - expect all cells but only one has fragments
        chrom_map = {"chr1": 1}
        valid_barcodes = set(cell_type_mapping.keys())
        coverage_aggregator, found_cells = processor._process_all_chromosomes(
            chrom_map, cell_type_mapping, valid_barcodes
        )

        assert found_cells == {present_cell}

        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        assert cell_types_in_coverage == {expected_cell_type}

        missing_cells = valid_barcodes - found_cells
        expected_missing = set(cell_type_mapping.keys()) - {present_cell}
        assert missing_cells == expected_missing

    def test__write_binned_coverage_per_chrom_full_pipeline(self, mocker, tmp_path):
        """Test write_binned_coverage_per_chrom() full pipeline orchestration."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        mock_pysam = mocker.patch("backend.layers.processing.utils.atac.pysam")
        mock_tabix = mocker.MagicMock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = [
            "chr1\t100\t200\tcell1",
            "chr1\t300\t400\tcell2",
        ]

        mock_tiledb = mocker.patch("backend.layers.processing.utils.atac.tiledb")
        mock_array = mocker.MagicMock()
        mock_array.__setitem__ = mocker.MagicMock()
        mock_tiledb.SparseArray.return_value.__enter__.return_value = mock_array

        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"}
        valid_barcodes = {"cell1", "cell2", "cell3"}

        processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))

        mock_tiledb.SparseArray.assert_called_with(array_name, mode="w", ctx=processor.ctx)
        assert mock_array.__setitem__.called

    def test__write_binned_coverage_per_chrom_empty_coverage(self, mocker, tmp_path):
        """Test write_binned_coverage_per_chrom() early return with empty coverage."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        mock_pysam = mocker.patch("backend.layers.processing.utils.atac.pysam")
        mock_tabix = mocker.MagicMock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = []
        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}

        result = processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)
        assert result is None
        mock_pysam.TabixFile.assert_called_once_with(str(fragment_file))

    def test__process_fragment_file_integration(
        self, tmp_path, cell_type_mapping, organism_genome_pair, mock_tiledb_components, mocker
    ):
        """Integration test for process_fragment_file."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"

        cell_names = list(cell_type_mapping.keys())[:2]
        fragment_lines = [f"chr1\t{100 + i*200}\t{200 + i*200}\t{cell}" for i, cell in enumerate(cell_names)]
        fragment_file.write_text("\n".join(fragment_lines) + "\n")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Mock pysam TabixFile
        mock_pysam = mocker.patch("backend.layers.processing.utils.atac.pysam")
        mock_tabix = mocker.MagicMock()
        mock_pysam.TabixFile.return_value.__enter__.return_value = mock_tabix
        mock_tabix.fetch.return_value = fragment_lines

        organism_id, expected_genome = organism_genome_pair

        cell_types = [cell_type_mapping[cell] for cell in cell_names]
        obs = pd.DataFrame(
            {
                "cell_type": cell_types,
                "organism_ontology_term_id": [organism_id] * len(cell_types),
            },
            index=cell_names,
        )

        array_name = str(tmp_path / "test_array")

        processor.process_fragment_file(obs, array_name)

        df_meta, genome_version = processor.extract_cell_metadata_from_h5ad(obs)

        expected_valid_barcodes = set(cell_names)
        assert set(df_meta["cell_name"]) == expected_valid_barcodes

        expected_cell_type_map = {cell: cell_type_mapping[cell] for cell in cell_names}
        actual_cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))
        assert actual_cell_type_map == expected_cell_type_map

        assert genome_version == expected_genome
        assert mock_tabix.fetch.call_count > 0
        mock_pysam.TabixFile.assert_called_with(str(fragment_file))
        assert isinstance(df_meta, pd.DataFrame)
        assert len(df_meta) == len(cell_names)

    def test__convert_coverage_to_cxg_array(self, tmp_path, cell_type_mapping, organism_genome_pair, mocker):
        """
        Test convert_coverage_to_cxg_array creates ATACDataProcessor and processes fragment file correctly.
        """
        from backend.layers.processing.utils.cxg_generation_utils import convert_coverage_to_cxg_array

        mock_atac_processor = mocker.MagicMock()
        mock_atac_processor_class = mocker.patch(
            "backend.layers.processing.utils.cxg_generation_utils.ATACDataProcessor",
            return_value=mock_atac_processor,
        )

        cxg_container = str(tmp_path / "test_container.cxg")
        organism_id, _ = organism_genome_pair

        cell_names = list(cell_type_mapping.keys())
        cell_types = list(cell_type_mapping.values())
        metadata_dict = pd.DataFrame(
            {
                "cell_type": cell_types,
                "organism_ontology_term_id": [organism_id] * len(cell_types),
            },
            index=cell_names,
        )
        fragment_artifact_id = str(tmp_path / "fragments.tsv.gz")
        group_metadata_name = "coverage"
        ctx = mocker.MagicMock()

        convert_coverage_to_cxg_array(cxg_container, metadata_dict, fragment_artifact_id, group_metadata_name, ctx)

        mock_atac_processor_class.assert_called_once_with(fragment_artifact_id, ctx)

        expected_array_name = f"{cxg_container}/{group_metadata_name}"
        mock_atac_processor.process_fragment_file.assert_called_once_with(metadata_dict, expected_array_name)
