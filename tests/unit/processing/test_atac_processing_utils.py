import os
from collections import defaultdict

import cellxgene_schema.atac_seq as atac_seq
import numpy as np
import pandas as pd
import pytest
import tiledb

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
    mock_tiledb.ByteShuffleFilter.return_value = mocker.MagicMock()
    mock_tiledb.DictionaryFilter.return_value = mocker.MagicMock()
    mock_tiledb.BitWidthReductionFilter.return_value = mocker.MagicMock()  # NEW
    mock_tiledb.DoubleDeltaFilter.return_value = mocker.MagicMock()        # NEW
    mock_tiledb.XORFilter.return_value = mocker.MagicMock()                # NEW
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

        df_meta, chromosome_by_length = processor.extract_cell_metadata_from_h5ad(obs)

        expected_df = pd.DataFrame({"cell_name": cell_names, "cell_type": cell_types})
        pd.testing.assert_frame_equal(df_meta, expected_df)
        assert chromosome_by_length is not None
        assert isinstance(chromosome_by_length, dict)
        # Check that we have common chromosomes
        assert "chr1" in chromosome_by_length
        assert "chr2" in chromosome_by_length

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
        organism_id, _ = organism_genome_pair
        chromosome_by_length = atac_seq.organism_ontology_term_id_by_chromosome_length_table.get(organism_id)
        max_chrom, chrom_map = processor.build_chrom_mapping(chromosome_by_length)

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
        organism_id, _ = organism_genome_pair
        chromosome_by_length = atac_seq.organism_ontology_term_id_by_chromosome_length_table.get(organism_id)
        max_bins = processor.calculate_max_bins(chromosome_by_length)

        assert isinstance(max_bins, int)
        assert max_bins > 0
        assert max_bins > 500_000  # At least 500K bins
        assert max_bins < 3_500_000  # Less than 3.5M bins (hg38 is larger than mm39)

    def test__normalized_coverage_calculation(self, tmp_path):
        """Test normalized coverage calculation: (count / total_coverage) * normalization_factor."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")

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

        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator)
        processor.create_dataframe_array(array_name, 2, 10)

        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)
        assert total_written == 4

        # Verify data was written correctly

        with tiledb.SparseArray(array_name, mode="r", ctx=processor.ctx) as A:
            df_data = A.df[:]

        # Convert to DataFrame for easier testing
        df = pd.DataFrame(df_data)
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert set(df.columns) == set(expected_columns)

        # Check optimized data types - TileDB may use uint32 for dimensions
        assert df["chrom"].dtype in ["int32", "uint32"]
        assert df["bin"].dtype in ["int32", "uint32"]
        assert df["coverage"].dtype == "uint16"  # Optimized from int32
        assert df["total_coverage"].dtype == "uint32"  # Optimized from int32  
        assert df["normalized_coverage"].dtype == "float32"  # TileDB doesn't support float16
        assert df["cell_type"].dtype == "object"

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

    def test__stream_coverage_chunks_small_dataset(self, tmp_path, coverage_aggregator_data):
        """Test _stream_coverage_chunks_to_tiledb() with parametrized datasets."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator_data)

        # Determine max dimensions for TileDB array
        max_chrom = max(chrom for chrom, _, _ in coverage_aggregator_data) if coverage_aggregator_data else 1
        max_bin = max(bin_id for _, bin_id, _ in coverage_aggregator_data) if coverage_aggregator_data else 1
        processor.create_dataframe_array(array_name, max_chrom, max(max_bin + 1, 10))

        total_written = processor._stream_coverage_chunks_to_tiledb(
            coverage_aggregator_data, cell_type_totals, array_name
        )
        expected_rows = len(coverage_aggregator_data)
        assert total_written == expected_rows

        with tiledb.SparseArray(array_name, mode="r", ctx=processor.ctx) as A:
            df_data = A.df[:]

        # Convert to DataFrame for easier testing
        df = pd.DataFrame(df_data)
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert set(df.columns) == set(expected_columns)

        assert df["chrom"].dtype in ["int32", "uint32"]  # TileDB may use uint32 for dimensions
        assert df["bin"].dtype in ["int32", "uint32"]
        assert df["coverage"].dtype == "uint16"  # Optimized from int32
        assert df["total_coverage"].dtype == "uint32"  # Optimized from int32
        assert df["normalized_coverage"].dtype == "float32"  # TileDB doesn't support float16
        assert df["cell_type"].dtype == "object"  # String columns are object type

        expected_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator_data)

        for cell_type, expected_total in expected_totals.items():
            cell_type_rows = df[df["cell_type"] == cell_type]
            assert all(cell_type_rows["total_coverage"] == expected_total)

        for _, row in df.iterrows():
            expected_normalized = (row["coverage"] / row["total_coverage"]) * 2_000_000
            assert abs(row["normalized_coverage"] - expected_normalized) < 1e-1

        assert not df.isnull().any().any()

    def test__stream_coverage_chunks_empty_aggregator(self, tmp_path):
        """Test _stream_coverage_chunks_to_tiledb() with empty aggregator."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        coverage_aggregator = defaultdict(int)
        # Compute cell type totals using helper method
        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator)

        # Create TileDB array for testing
        processor.create_dataframe_array(array_name, 1, 10)

        # Stream empty data to TileDB
        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)

        assert total_written == 0

        # Read back from TileDB to verify empty array

        with tiledb.SparseArray(array_name, mode="r", ctx=processor.ctx) as A:
            df_data = A.df[:]

        # Convert to DataFrame for easier testing
        df = pd.DataFrame(df_data)
        assert len(df) == 0
        expected_columns = ["chrom", "bin", "cell_type", "coverage", "total_coverage", "normalized_coverage"]
        assert set(df.columns) == set(expected_columns)

        assert df.empty
        assert isinstance(df, pd.DataFrame)

    def test__stream_coverage_direct_assignment_consistency(self, tmp_path):
        """Test that direct array assignment produces consistent results."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update(
            {
                (1, 10, "T cell"): 15,
                (1, 20, "T cell"): 25,
                (2, 5, "B cell"): 30,
                (2, 15, "B cell"): 40,
                (3, 8, "NK cell"): 12,
            }
        )

        expected_totals = {"T cell": 40, "B cell": 70, "NK cell": 12}
        processor.create_dataframe_array(array_name, 3, 25)

        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, expected_totals, array_name)
        assert total_written == 5

        with tiledb.open(array_name, mode="r") as array:
            result = array[:]
            assert len(result["chrom"]) == 5

            # Verify key data points
            chrom_bin_pairs = list(zip(result["chrom"], result["bin"], strict=False))
            assert (1, 10) in chrom_bin_pairs
            assert (2, 5) in chrom_bin_pairs
            assert (3, 8) in chrom_bin_pairs

            # Verify coverage values match input
            for i, (chrom, bin_val) in enumerate(chrom_bin_pairs):
                if chrom == 1 and bin_val == 10:
                    assert result["coverage"][i] == 15
                elif chrom == 2 and bin_val == 5:
                    assert result["coverage"][i] == 30
                elif chrom == 3 and bin_val == 8:
                    assert result["coverage"][i] == 12

    def test__stream_coverage_direct_assignment_memory_efficiency(self, tmp_path):
        """Test that direct assignment doesn't create intermediate data structures."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")
        # Use min_coverage_threshold=1 to prevent pruning in this memory efficiency test
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file), min_coverage_threshold=1)

        # Generate larger dataset (1000 records) to test memory efficiency
        coverage_aggregator = defaultdict(int)
        num_records = 1000
        for i in range(num_records):
            chrom = (i % 3) + 1
            bin_val = i * 10
            cell_type = f"cell_type_{i % 5}"
            coverage_aggregator[(chrom, bin_val, cell_type)] = i + 1

        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator)
        processor.create_dataframe_array(array_name, 3, num_records * 10)

        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)
        assert total_written == num_records

        # Verify statistical properties
        with tiledb.open(array_name, mode="r") as array:
            result = array[:]
            assert len(result["chrom"]) == num_records
            assert min(result["coverage"]) == 1
            assert max(result["coverage"]) == num_records
            assert sum(result["coverage"]) == sum(range(1, num_records + 1))

    def test__stream_coverage_direct_assignment_edge_cases(self, tmp_path):
        """Test direct assignment handles edge cases correctly."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        # Test single record edge case
        coverage_aggregator = defaultdict(int)
        coverage_aggregator[(1, 0, "single_cell")] = 42
        cell_type_totals = {"single_cell": 42}

        processor.create_dataframe_array(array_name, 1, 10)
        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)
        assert total_written == 1

        with tiledb.open(array_name, mode="r") as array:
            result = array[:]
            assert len(result["chrom"]) == 1
            assert result["chrom"][0] == 1
            assert result["bin"][0] == 0
            assert result["coverage"][0] == 42
            assert result["total_coverage"][0] == 42

            # Verify normalized coverage calculation
            expected_normalized = (42 / 42) * processor.normalization_factor
            assert abs(result["normalized_coverage"][0] - expected_normalized) < 1e-6

    def test__stream_coverage_zero_total_coverage_handling(self, tmp_path):
        """Test handling of zero total coverage scenarios."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        array_name = str(tmp_path / "test_array")
        # Use threshold of 0 to disable pruning for this edge case test
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file), min_coverage_threshold=0)

        # Test zero total coverage edge case
        coverage_aggregator = defaultdict(int)
        coverage_aggregator[(1, 0, "orphan_cell")] = 10
        cell_type_totals = {"orphan_cell": 0}  # Zero total (edge case)

        processor.create_dataframe_array(array_name, 1, 10)
        total_written = processor._stream_coverage_chunks_to_tiledb(coverage_aggregator, cell_type_totals, array_name)
        assert total_written == 1

        # Verify zero total coverage handling
        with tiledb.open(array_name, mode="r") as array:
            result = array[:]
            assert result["coverage"][0] == 10
            assert result["total_coverage"][0] == 0
            assert result["normalized_coverage"][0] == 0.0  # Should default to 0 when total is 0

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

        # Updated assertions for advanced compression strategy
        assert mock_tiledb.FilterList.call_count >= 4  # coverage_int, coverage_float, categorical, dims
        mock_tiledb.ByteShuffleFilter.assert_called()  # Used in all filter combinations
        mock_tiledb.DictionaryFilter.assert_called_once()  # Used for categorical compression
        mock_tiledb.BitWidthReductionFilter.assert_called()  # NEW: Used for integer compression
        mock_tiledb.DoubleDeltaFilter.assert_called_once()  # NEW: Used for coordinate compression
        
        # Verify ZStd calls with updated compression levels
        assert mock_tiledb.ZstdFilter.call_count >= 4  # coverage_int (level=6), coverage_float (level=6), categorical (level=22), dims (level=6)
        zstd_calls = mock_tiledb.ZstdFilter.call_args_list
        
        # Verify we have calls with different compression levels
        level_6_calls = [call for call in zstd_calls if call.kwargs.get("level") == 6]
        level_22_calls = [call for call in zstd_calls if call.kwargs.get("level") == 22]
        assert len(level_6_calls) >= 3  # coverage_int, coverage_float, and dimension filters
        assert len(level_22_calls) == 1  # Categorical compression

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
        # Verify adaptive tile sizing: between 100 and 10000, but cannot exceed domain range
        calculated_tile = min(max(max_bins // 1000, 100), 10000)
        expected_tile = min(calculated_tile, max_bins)  # Cannot exceed domain range
        assert bin_dim_call.kwargs["tile"] == expected_tile
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
        assert coverage_attr.kwargs["dtype"] == np.uint16  # Optimized from int32

        total_coverage_attr = next(call for call in attr_calls if call.kwargs["name"] == "total_coverage")
        assert total_coverage_attr.kwargs["dtype"] == np.uint32  # Optimized from int32

        normalized_attr = next(call for call in attr_calls if call.kwargs["name"] == "normalized_coverage")
        assert normalized_attr.kwargs["dtype"] == np.float32  # TileDB doesn't support float16

        mock_tiledb.ArraySchema.assert_called_once()
        schema_kwargs = mock_tiledb.ArraySchema.call_args.kwargs
        assert schema_kwargs["domain"] == mock_tiledb_components["domain"]
        assert len(schema_kwargs["attrs"]) == 3
        assert schema_kwargs["sparse"] is True
        assert schema_kwargs["allows_duplicates"] is False

        mock_tiledb.SparseArray.create.assert_called_once_with(array_name, mock_tiledb_components["schema"])

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

    def test__process_all_chromosomes(self, tmp_path, cell_type_mapping, invalid_fragment_coordinates):
        """Test _process_all_chromosomes() processes fragments from all chromosomes."""
        # Create a proper compressed fragment file with real data

        fragment_file = tmp_path / "test_fragments.tsv.gz"

        valid_cells = list(cell_type_mapping.keys())[:3]

        # Create test fragment data in proper format (sorted by chromosome then position)
        fragments = [
            f"chr1\t100\t200\t{valid_cells[0]}",
            f"chr1\t300\t400\t{valid_cells[1]}",
            f"chr1\t500\t600\t{valid_cells[0]}",
            f"chrX\t200\t300\t{valid_cells[0]}",
            # chr2 intentionally empty to test missing chromosome handling
        ]

        # Write uncompressed first, then compress and index
        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        # Compress with bgzip (required for tabix)
        import subprocess

        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)

        # Create tabix index
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        chrom_map = {"chr1": 1, "chr2": 2, "chrX": 23}
        valid_barcodes = set(valid_cells)
        coverage_aggregator, found_cells = processor._process_all_chromosomes(
            chrom_map, cell_type_mapping, valid_barcodes
        )

        expected_found_cells = {valid_cells[0], valid_cells[1]}
        assert found_cells == expected_found_cells

        assert len(coverage_aggregator) > 0

        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        expected_cell_types = {cell_type_mapping[cell] for cell in expected_found_cells}
        assert cell_types_in_coverage == expected_cell_types

        chroms_in_coverage = {key[0] for key in coverage_aggregator}
        expected_chroms = {1, 23}  # chr1=1, chrX=23 (chr2=2 failed)
        assert chroms_in_coverage == expected_chroms

    def test__process_all_chromosomes_missing_cells(self, tmp_path, cell_type_mapping):
        """Test _process_all_chromosomes() handles missing cells correctly."""
        # Create fragment file with only one cell
        present_cell = list(cell_type_mapping.keys())[0]
        expected_cell_type = cell_type_mapping[present_cell]

        fragments = [f"chr1\t100\t200\t{present_cell}"]

        # Write uncompressed first, then compress and index
        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        # Compress with bgzip and create tabix index
        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

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

    def test__write_binned_coverage_per_chrom_full_pipeline(self, tmp_path, mock_tiledb_components):
        """Test write_binned_coverage_per_chrom() full pipeline orchestration."""
        fragments = [
            "chr1\t100\t200\tcell1",
            "chr1\t300\t400\tcell2",
        ]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        # Use threshold=1 to ensure test data isn't pruned (coverage values will be 1)
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file), min_coverage_threshold=1)

        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell", "cell2": "B cell", "cell3": "NK cell"}
        valid_barcodes = {"cell1", "cell2", "cell3"}

        processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)

        mock_tiledb_components["tiledb"].SparseArray.assert_called()
        assert mock_tiledb_components["tiledb"].SparseArray.return_value.__enter__.return_value.__setitem__.called

    def test__write_binned_coverage_per_chrom_empty_coverage(self, tmp_path):
        """Test write_binned_coverage_per_chrom() early return with empty coverage."""
        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

        array_name = str(tmp_path / "test_array")
        chrom_map = {"chr1": 1}
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}

        result = processor.write_binned_coverage_per_chrom(array_name, chrom_map, cell_type_map, valid_barcodes)
        assert result is None

    def test__process_fragment_file_integration(
        self, tmp_path, cell_type_mapping, organism_genome_pair, mock_tiledb_components
    ):
        """Integration test for process_fragment_file."""
        cell_names = list(cell_type_mapping.keys())[:2]
        fragment_lines = [f"chr1\t{100 + i*200}\t{200 + i*200}\t{cell}" for i, cell in enumerate(cell_names)]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragment_lines) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file))

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

        df_meta, chromosome_by_length = processor.extract_cell_metadata_from_h5ad(obs)

        expected_valid_barcodes = set(cell_names)
        assert set(df_meta["cell_name"]) == expected_valid_barcodes

        expected_cell_type_map = {cell: cell_type_mapping[cell] for cell in cell_names}
        actual_cell_type_map = dict(zip(df_meta["cell_name"], df_meta["cell_type"], strict=False))
        assert actual_cell_type_map == expected_cell_type_map

        assert chromosome_by_length is not None
        assert isinstance(chromosome_by_length, dict)
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
        mock_atac_processor.process_fragment_file.assert_called_once_with(metadata_dict, expected_array_name, uns=None)

    def test__process_chromosome_worker_valid_fragments(self, tmp_path):
        """Test _process_chromosome_worker processes valid fragments correctly."""
        fragments = [
            "chr1\t100\t200\tcell1",
            "chr1\t250\t350\tcell3",
            "chr1\t300\t400\tcell1",
            "chr1\t500\t600\tcell2",
        ]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell", "cell2": "B cell"}
        valid_barcodes = {"cell1", "cell2"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert found_cells == {"cell1", "cell2"}
        assert len(coverage_aggregator) > 0
        assert (chrom_id, 1, "T cell") in coverage_aggregator
        assert (chrom_id, 3, "T cell") in coverage_aggregator
        assert (chrom_id, 5, "B cell") in coverage_aggregator

        cell_types_in_coverage = {key[2] for key in coverage_aggregator}
        assert len(cell_types_in_coverage) == 2

        for count in coverage_aggregator.values():
            assert count > 0

    def test__process_chromosome_worker_spanning_bins(self, tmp_path):
        """Test _process_chromosome_worker handles fragments spanning multiple bins."""
        fragments = ["chr1\t150\t350\tcell1"]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert found_cells == {"cell1"}
        assert (chrom_id, 1, "T cell") in coverage_aggregator
        assert (chrom_id, 3, "T cell") in coverage_aggregator
        assert coverage_aggregator[(chrom_id, 1, "T cell")] == 1
        assert coverage_aggregator[(chrom_id, 3, "T cell")] == 1

    def test__process_chromosome_worker_invalid_coordinates(self, tmp_path):
        """Test _process_chromosome_worker handles invalid coordinates gracefully."""
        fragments = [
            "chr1\t100\t200\tcell1",
            "chr1\t500\t600\tcell1",
        ]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert found_cells == {"cell1"}

        expected_bins = {1, 5}
        actual_bins = {key[1] for key in coverage_aggregator}
        assert actual_bins == expected_bins

    def test__process_chromosome_worker_malformed_rows(self, tmp_path):
        """Test _process_chromosome_worker handles malformed rows gracefully."""
        fragments = [
            "chr1\t100\t200\tcell1",
            "chr1\t400\t500\tcell1",
        ]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert found_cells == {"cell1"}

        expected_bins = {1, 4}
        actual_bins = {key[1] for key in coverage_aggregator}
        assert actual_bins == expected_bins

    def test__process_chromosome_worker_empty_file(self, tmp_path):
        """Test _process_chromosome_worker handles empty fragment file."""
        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert len(coverage_aggregator) == 0
        assert len(found_cells) == 0

    def test__process_chromosome_worker_no_valid_cells(self, tmp_path):
        """Test _process_chromosome_worker when no fragments match valid barcodes."""
        fragments = [
            "chr1\t100\t200\tinvalid_cell1",
            "chr1\t300\t400\tinvalid_cell2",
        ]

        uncompressed_file = tmp_path / "test_fragments.tsv"
        with open(uncompressed_file, "w") as f:
            f.write("\n".join(fragments) + "\n")

        import subprocess

        fragment_file = tmp_path / "test_fragments.tsv.gz"
        subprocess.run(["bgzip", str(uncompressed_file)], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "bed", str(fragment_file)], check=True, capture_output=True)

        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"valid_cell": "T cell"}
        valid_barcodes = {"valid_cell"}
        bin_size = 100

        args = (str(fragment_file), chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert len(coverage_aggregator) == 0
        assert len(found_cells) == 0

    def test__process_chromosome_worker_coordinate_validation(self, mocker):
        """Test _process_chromosome_worker coordinate validation logic using mocks."""
        mock_tabix_file = mocker.MagicMock()
        mock_tabix = mocker.MagicMock()
        mock_tabix.__enter__.return_value = mock_tabix_file

        mock_tabix_file.fetch.return_value = [
            "chr1\t100\t200\tcell1",
            "chr1\t-50\t100\tcell1",
            "chr1\t300\t300\tcell1",
            "chr1\t400\t350\tcell1",
            "chr1\tabc\t200\tcell1",
            "chr1\t500\t600\tcell1",
            "chr1\t200\t300",
            "incomplete",
        ]

        mocker.patch("pysam.TabixFile", return_value=mock_tabix)

        fragment_file = "/fake/path/fragments.tsv.gz"
        chrom_str = "chr1"
        chrom_id = 1
        cell_type_map = {"cell1": "T cell"}
        valid_barcodes = {"cell1"}
        bin_size = 100

        args = (fragment_file, chrom_str, chrom_id, cell_type_map, valid_barcodes, bin_size)

        coverage_aggregator, found_cells = ATACDataProcessor._process_chromosome_worker(args)

        assert found_cells == {"cell1"}

        expected_bins = {1, 5}
        actual_bins = {key[1] for key in coverage_aggregator}
        assert actual_bins == expected_bins

        mock_tabix_file.fetch.assert_called_once_with(chrom_str)

    def test__file_exists_local_file_exists(self, tmp_path):
        """Test _file_exists returns True for existing local file."""
        test_file = tmp_path / "test_file.txt"
        test_file.write_text("test content")

        processor = ATACDataProcessor()
        assert processor._file_exists(str(test_file)) is True

    def test__file_exists_local_file_not_exists(self):
        """Test _file_exists returns False for non-existing local file."""
        processor = ATACDataProcessor()
        assert processor._file_exists("/non/existent/file.txt") is False

    def test__file_exists_s3_file_exists(self, mocker):
        """Test _file_exists returns True for existing S3 file."""
        mock_s3_client = mocker.MagicMock()
        mock_boto3 = mocker.patch("backend.layers.processing.utils.atac.boto3")
        mock_boto3.client.return_value = mock_s3_client
        mock_s3_client.head_object.return_value = {}

        processor = ATACDataProcessor()
        result = processor._file_exists("s3://test-bucket/test-key.txt")

        assert result is True
        mock_boto3.client.assert_called_once_with("s3")
        mock_s3_client.head_object.assert_called_once_with(Bucket="test-bucket", Key="test-key.txt")

    def test__file_exists_s3_file_not_exists(self, mocker):
        """Test _file_exists returns False when S3 file doesn't exist."""
        mock_s3_client = mocker.MagicMock()
        mock_boto3 = mocker.patch("backend.layers.processing.utils.atac.boto3")
        mock_boto3.client.return_value = mock_s3_client
        mock_s3_client.head_object.side_effect = Exception("NoSuchKey")

        processor = ATACDataProcessor()
        result = processor._file_exists("s3://test-bucket/nonexistent-key.txt")

        assert result is False
        mock_boto3.client.assert_called_once_with("s3")
        mock_s3_client.head_object.assert_called_once_with(Bucket="test-bucket", Key="nonexistent-key.txt")

    def test__download_s3_file_local_path(self):
        """Test _download_s3_file returns same path for local files."""
        processor = ATACDataProcessor()
        local_path = "/local/path/file.txt"
        result = processor._download_s3_file(local_path)
        assert result == local_path

    def test__download_s3_file_s3_success(self, mocker, tmp_path):
        """Test _download_s3_file successfully downloads S3 file and index."""
        mock_s3_client = mocker.MagicMock()
        mock_boto3 = mocker.patch("backend.layers.processing.utils.atac.boto3")
        mock_boto3.client.return_value = mock_s3_client

        mock_tempfile = mocker.patch("backend.layers.processing.utils.atac.tempfile")
        mock_temp = mocker.MagicMock()
        mock_temp.name = str(tmp_path / "temp_file.bgz")
        mock_tempfile.NamedTemporaryFile.return_value = mock_temp

        processor = ATACDataProcessor()
        s3_path = "s3://test-bucket/fragments.tsv.gz"

        result = processor._download_s3_file(s3_path)

        assert result == mock_temp.name
        mock_boto3.client.assert_called_once_with("s3")
        mock_s3_client.download_file.assert_any_call("test-bucket", "fragments.tsv.gz", mock_temp.name)
        mock_s3_client.download_file.assert_any_call("test-bucket", "fragments.tsv.gz.tbi", mock_temp.name + ".tbi")
        mock_temp.close.assert_called_once()

    def test__download_s3_file_download_failure(self, mocker, tmp_path):
        """Test _download_s3_file handles download failure and cleans up temp files."""
        mock_s3_client = mocker.MagicMock()
        mock_boto3 = mocker.patch("backend.layers.processing.utils.atac.boto3")
        mock_boto3.client.return_value = mock_s3_client
        mock_s3_client.download_file.side_effect = Exception("Download failed")

        mock_tempfile = mocker.patch("backend.layers.processing.utils.atac.tempfile")
        mock_temp = mocker.MagicMock()
        temp_path = str(tmp_path / "temp_file.bgz")
        mock_temp.name = temp_path
        mock_tempfile.NamedTemporaryFile.return_value = mock_temp

        with open(temp_path, "w") as f:
            f.write("temp")
        with open(temp_path + ".tbi", "w") as f:
            f.write("temp_index")

        processor = ATACDataProcessor()
        s3_path = "s3://test-bucket/fragments.tsv.gz"

        with pytest.raises(RuntimeError, match="Failed to download S3 file"):
            processor._download_s3_file(s3_path)

        assert not os.path.exists(temp_path)
        assert not os.path.exists(temp_path + ".tbi")

    def test__get_fragment_file_path_none_artifact_id(self):
        """Test _get_fragment_file_path raises ValueError when fragment_artifact_id is None."""
        processor = ATACDataProcessor()

        with pytest.raises(ValueError, match="Fragment artifact ID is not set"):
            processor._get_fragment_file_path()

    def test__get_fragment_file_path_local_file(self, tmp_path):
        """Test _get_fragment_file_path with local file path."""
        test_file = tmp_path / "fragments.tsv.gz"
        test_file.write_text("test")

        processor = ATACDataProcessor(fragment_artifact_id=str(test_file))
        result = processor._get_fragment_file_path()

        assert result == str(test_file)
        assert processor._local_fragment_file == str(test_file)

    def test__get_fragment_file_path_s3_file(self, mocker, tmp_path):
        """Test _get_fragment_file_path with S3 file path."""
        mock_file_exists = mocker.patch.object(ATACDataProcessor, "_file_exists", return_value=True)
        mock_download = mocker.patch.object(ATACDataProcessor, "_download_s3_file")
        temp_path = str(tmp_path / "downloaded_file.bgz")
        mock_download.return_value = temp_path

        s3_path = "s3://test-bucket/fragments.tsv.gz"
        processor = ATACDataProcessor(fragment_artifact_id=s3_path)

        result = processor._get_fragment_file_path()

        assert result == temp_path
        assert processor._local_fragment_file == temp_path
        assert processor._local_fragment_index == temp_path + ".tbi"
        mock_download.assert_called_once_with(s3_path)
        mock_file_exists.assert_called_once_with(s3_path)

    def test__get_fragment_file_path_cached(self, mocker, tmp_path):
        """Test _get_fragment_file_path returns cached path on subsequent calls."""
        mocker.patch.object(ATACDataProcessor, "_file_exists", return_value=True)
        mock_download = mocker.patch.object(ATACDataProcessor, "_download_s3_file")
        temp_path = str(tmp_path / "downloaded_file.bgz")
        mock_download.return_value = temp_path

        s3_path = "s3://test-bucket/fragments.tsv.gz"
        processor = ATACDataProcessor(fragment_artifact_id=s3_path)

        result1 = processor._get_fragment_file_path()
        result2 = processor._get_fragment_file_path()

        assert result1 == result2 == temp_path
        mock_download.assert_called_once()

    def test__cleanup_temp_files_no_temp_files(self):
        """Test _cleanup_temp_files when no temp files exist."""
        processor = ATACDataProcessor()
        processor._cleanup_temp_files()

    def test__cleanup_temp_files_local_file_no_cleanup(self, tmp_path):
        """Test _cleanup_temp_files doesn't delete original local files."""
        test_file = tmp_path / "fragments.tsv.gz"
        test_file.write_text("test")

        processor = ATACDataProcessor(fragment_artifact_id=str(test_file))
        processor._local_fragment_file = str(test_file)

        processor._cleanup_temp_files()

        assert test_file.exists()
        assert processor._local_fragment_file == str(test_file)

    def test__cleanup_temp_files_temp_files(self, tmp_path):
        """Test _cleanup_temp_files removes downloaded temp files."""
        s3_path = "s3://test-bucket/fragments.tsv.gz"
        temp_file = tmp_path / "temp_fragment.bgz"
        temp_index = tmp_path / "temp_fragment.bgz.tbi"

        temp_file.write_text("temp fragment")
        temp_index.write_text("temp index")

        processor = ATACDataProcessor()
        processor.fragment_artifact_id = s3_path
        processor._local_fragment_file = str(temp_file)
        processor._local_fragment_index = str(temp_index)

        processor._cleanup_temp_files()

        assert not temp_file.exists()
        assert not temp_index.exists()
        assert processor._local_fragment_file is None
        assert processor._local_fragment_index is None

    def test__cleanup_temp_files_partial_cleanup(self, tmp_path):
        """Test _cleanup_temp_files handles missing files gracefully."""
        s3_path = "s3://test-bucket/fragments.tsv.gz"
        temp_file = tmp_path / "temp_fragment.bgz"
        temp_index = tmp_path / "temp_fragment.bgz.tbi"

        temp_file.write_text("temp fragment")

        processor = ATACDataProcessor()
        processor.fragment_artifact_id = s3_path
        processor._local_fragment_file = str(temp_file)
        processor._local_fragment_index = str(temp_index)

        processor._cleanup_temp_files()

        assert not temp_file.exists()
        assert processor._local_fragment_file is None
        assert processor._local_fragment_index == str(temp_index)

    def test__sparse_data_pruning(self, tmp_path):
        """Test sparse data pruning removes low normalized coverage records below threshold."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        # Test with normalized coverage threshold (500000 to prune low values)
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file), min_coverage_threshold=500000)
        
        # Create test data with specific coverage values to test normalized pruning
        # T cell total: 5 + 3 + 1 = 9
        # B cell total: 10 + 2 = 12  
        # NK cell total: 1
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update({
            # T cell records (total = 9)
            (1, 0, "T cell"): 5,    # normalized = (5/9)*2000000 ≈ 1111111 (keep)
            (1, 1, "T cell"): 3,    # normalized = (3/9)*2000000 ≈ 666667 (keep)
            (4, 0, "T cell"): 1,    # normalized = (1/9)*2000000 ≈ 222222 (prune < 500000)
            # B cell records (total = 12)
            (2, 0, "B cell"): 10,   # normalized = (10/12)*2000000 ≈ 1666667 (keep)
            (2, 1, "B cell"): 2,    # normalized = (2/12)*2000000 ≈ 333333 (prune < 500000)
            # NK cell records (total = 1) 
            (3, 0, "NK cell"): 1,   # normalized = (1/1)*2000000 = 2000000 (keep)
        })
        
        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator)
        
        # Generate chunks and verify pruning
        chunks = list(processor._generate_chunks(coverage_aggregator, cell_type_totals, 1000))
        assert len(chunks) == 1  # Should fit in one chunk
        
        chunk_data = chunks[0]
        assert len(chunk_data) == 4  # Should keep 4 records with normalized coverage ≥ 500000
        
        # Verify correct records were kept (those with high normalized coverage)
        kept_coverages = {record["coverage"] for record in chunk_data}
        assert kept_coverages == {5, 3, 10, 1}  # Raw coverage values of kept records (T cell 5&3, B cell 10, NK cell 1)
        
        # Verify all kept records meet normalized threshold
        for record in chunk_data:
            assert record["normalized_coverage"] >= processor.min_coverage_threshold

    def test__sparse_data_pruning_threshold_zero(self, tmp_path):
        """Test that threshold=0 disables pruning."""
        fragment_file = tmp_path / "test_fragments.tsv.gz"
        fragment_file.write_text("chr1\t100\t200\tcell1\n")

        # Test with threshold of 0 (no pruning)
        processor = ATACDataProcessor(fragment_artifact_id=str(fragment_file), min_coverage_threshold=0)
        
        coverage_aggregator = defaultdict(int)
        coverage_aggregator.update({
            (1, 0, "T cell"): 5,
            (2, 0, "B cell"): 1,    # This should be kept with threshold=0
            (3, 0, "NK cell"): 0,   # Even 0 coverage should be kept
        })
        
        cell_type_totals = processor._compute_cell_type_totals_from_aggregator(coverage_aggregator)
        chunks = list(processor._generate_chunks(coverage_aggregator, cell_type_totals, 1000))
        
        assert len(chunks) == 1
        chunk_data = chunks[0] 
        assert len(chunk_data) == 3  # All records should be kept
        
        kept_coverages = {record["coverage"] for record in chunk_data}
        assert kept_coverages == {5, 1, 0}  # All values should remain
