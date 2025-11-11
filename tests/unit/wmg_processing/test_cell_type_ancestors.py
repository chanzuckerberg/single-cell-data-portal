import json
import unittest
from unittest.mock import patch

from backend.common.census_cube.data.snapshot import CELL_TYPE_ANCESTORS_FILENAME
from backend.common.census_cube.utils import ancestors, children, descendants
from backend.wmg.pipeline.cell_type_ancestors import create_cell_type_ancestors
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ANCESTORS_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.utils import load_pipeline_state, write_pipeline_state
from tests.test_utils import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_tmpdir


class CellTypeAncestorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_cube_dir = load_realistic_test_snapshot_tmpdir("realistic-test-snapshot")
        with open(f"{cls.temp_cube_dir.name}/{CELL_TYPE_ANCESTORS_FILENAME}") as f:
            cls.expected_cell_type_ancestors = json.load(f)

    def setUp(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state[EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG] = True
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    def tearDown(self):
        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        pipeline_state.pop(CELL_TYPE_ANCESTORS_CREATED_FLAG, None)
        write_pipeline_state(pipeline_state, self.temp_cube_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_cube_dir.cleanup()

    def test_cell_type_ancestors(self):
        create_cell_type_ancestors(self.temp_cube_dir.name)

        with open(f"{self.temp_cube_dir.name}/{CELL_TYPE_ANCESTORS_FILENAME}") as f:
            cell_type_ancestors = json.load(f)

        pipeline_state = load_pipeline_state(self.temp_cube_dir.name)
        self.assertTrue(pipeline_state.get(CELL_TYPE_ANCESTORS_CREATED_FLAG))
        self.assertTrue(compare_dicts(cell_type_ancestors, self.expected_cell_type_ancestors))


class CellTypeUtilsTests(unittest.TestCase):
    """Unit tests for census_cube.utils functions."""

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_ancestors_handles_key_error(self, mock_ontology_parser):
        """Test that ancestors handles KeyError for missing cell types in ontology."""
        # Simulate the error that occurred in production with CL:4047053
        mock_ontology_parser.get_term_ancestors.side_effect = KeyError("CL:4047053")

        # Clear cache to ensure clean test
        ancestors.cache_clear()

        result = ancestors("CL:4047053")

        # Should return the cell type itself as a fallback instead of crashing
        self.assertEqual(result, ["CL:4047053"])

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_ancestors_handles_value_error(self, mock_ontology_parser):
        """Test that ancestors handles ValueError gracefully."""
        mock_ontology_parser.get_term_ancestors.side_effect = ValueError("Invalid term")

        # Clear cache to ensure clean test
        ancestors.cache_clear()

        result = ancestors("CL:INVALID")

        # Should return the cell type itself as a fallback
        self.assertEqual(result, ["CL:INVALID"])

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_descendants_handles_key_error(self, mock_ontology_parser):
        """Test that descendants handles KeyError for missing cell types in ontology."""
        mock_ontology_parser.get_term_descendants.side_effect = KeyError("CL:4047053")

        # Clear cache to ensure clean test
        descendants.cache_clear()

        result = descendants("CL:4047053")

        # Should return the cell type itself as a fallback instead of crashing
        self.assertEqual(result, ["CL:4047053"])

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_descendants_handles_value_error(self, mock_ontology_parser):
        """Test that descendants handles ValueError gracefully."""
        mock_ontology_parser.get_term_descendants.side_effect = ValueError("Invalid term")

        # Clear cache to ensure clean test
        descendants.cache_clear()

        result = descendants("CL:INVALID")

        # Should return the cell type itself as a fallback
        self.assertEqual(result, ["CL:INVALID"])

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_children_handles_key_error(self, mock_ontology_parser):
        """Test that children handles KeyError for missing cell types in ontology."""
        mock_ontology_parser.get_term_children.side_effect = KeyError("CL:4052026")

        # Clear cache to ensure clean test
        children.cache_clear()

        result = children("CL:4052026")

        # Should return an empty list as a fallback instead of crashing
        self.assertEqual(result, [])

    @patch("backend.common.census_cube.utils.ontology_parser")
    def test_children_handles_value_error(self, mock_ontology_parser):
        """Test that children handles ValueError gracefully."""
        mock_ontology_parser.get_term_children.side_effect = ValueError("Invalid term")

        # Clear cache to ensure clean test
        children.cache_clear()

        result = children("CL:INVALID")

        # Should return an empty list as a fallback
        self.assertEqual(result, [])
