import unittest
from unittest.mock import patch

from backend.common.census_cube.data.ontology_labels import (
    gene_term_label,
    is_ontology_term_deprecated,
    ontology_term_description,
    ontology_term_label,
    ontology_term_synonyms,
)


class OntologyLabelTests(unittest.TestCase):
    def test__ontology_term_label(self):
        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        # This test is fragile, in that changes to our ontologies may break these these tests in the future, but that
        # should happen infrequently
        test_cases = {
            "EFO:0000001": "experimental factor",
            "HANCESTRO:0003": "country",
            "CL:0000000": "cell",
            "HsapDv:0000000": "human life cycle stage",
            "PATO:0000001": "quality",
            "NCBITaxon:10001": "Sciurus",
            "MmusDv:0000000": "mouse life cycle stage",
            "MONDO:0000001": "disease",
            "UBERON:0002539": "pharyngeal arch",
        }
        for ontology_term_id, expected_ontology_term_label in test_cases.items():
            with self.subTest(ontology_term_id):
                self.assertEqual(ontology_term_label(ontology_term_id), expected_ontology_term_label)

    def test__gene_label(self):
        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        # This test is fragile, in that changes to our ontologies may break these these tests in the future, but that
        # should happen infrequently
        test_cases = {
            "ERCC-00002": "ERCC-00002 (spike-in control)",
            "ENSG00000223972": "DDX11L1",
            "ENSMUSG00000102693": "4933401J01Rik",
        }
        for gene_id, expected_gene_label in test_cases.items():
            with self.subTest(gene_id):
                self.assertEqual(gene_term_label(gene_id), expected_gene_label)

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_label_handles_key_error(self, mock_ontology_parser):
        """Test that ontology_term_label handles KeyError for missing terms in ontology."""
        # Simulate the error that occurred in production with CL:4033085
        mock_ontology_parser.get_term_label.side_effect = KeyError("CL:4033085")

        result = ontology_term_label("CL:4033085")

        # Should return the term ID itself as a fallback instead of crashing
        self.assertEqual(result, "CL:4033085")

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_label_handles_value_error(self, mock_ontology_parser):
        """Test that ontology_term_label handles ValueError gracefully."""
        mock_ontology_parser.get_term_label.side_effect = ValueError("Invalid term")

        result = ontology_term_label("CL:INVALID")

        # Should return the term ID itself as a fallback
        self.assertEqual(result, "CL:INVALID")

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_label_success_path(self, mock_ontology_parser):
        """Test that ontology_term_label returns the correct label when term exists."""
        mock_ontology_parser.get_term_label.return_value = "native cell"

        result = ontology_term_label("CL:0000003")

        self.assertEqual(result, "native cell")
        mock_ontology_parser.get_term_label.assert_called_once_with("CL:0000003")

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_is_ontology_term_deprecated_handles_key_error(self, mock_ontology_parser):
        """Test that is_ontology_term_deprecated handles KeyError for missing terms."""
        mock_ontology_parser.is_term_deprecated.side_effect = KeyError("CL:4033085")

        result = is_ontology_term_deprecated("CL:4033085")

        # Should return False for unknown terms (assume not deprecated)
        self.assertEqual(result, False)

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_is_ontology_term_deprecated_success_path(self, mock_ontology_parser):
        """Test that is_ontology_term_deprecated returns the correct value when term exists."""
        mock_ontology_parser.is_term_deprecated.return_value = True

        result = is_ontology_term_deprecated("CL:0000001")

        self.assertEqual(result, True)
        mock_ontology_parser.is_term_deprecated.assert_called_once_with("CL:0000001")

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_description_handles_key_error(self, mock_ontology_parser):
        """Test that ontology_term_description handles KeyError for missing terms."""
        mock_ontology_parser.get_term_description.side_effect = KeyError("CL:4033085")

        result = ontology_term_description("CL:4033085")

        # Should return None for unknown terms (no description available)
        self.assertIsNone(result)

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_description_success_path(self, mock_ontology_parser):
        """Test that ontology_term_description returns the correct description when term exists."""
        mock_ontology_parser.get_term_description.return_value = "A cell description"

        result = ontology_term_description("CL:0000003")

        self.assertEqual(result, "A cell description")
        mock_ontology_parser.get_term_description.assert_called_once_with("CL:0000003")

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_synonyms_handles_key_error(self, mock_ontology_parser):
        """Test that ontology_term_synonyms handles KeyError for missing terms."""
        mock_ontology_parser.get_term_synonyms.side_effect = KeyError("CL:4033085")

        result = ontology_term_synonyms("CL:4033085")

        # Should return empty list for unknown terms (no synonyms available)
        self.assertEqual(result, [])

    @patch("backend.common.census_cube.data.ontology_labels.ontology_parser")
    def test_ontology_term_synonyms_success_path(self, mock_ontology_parser):
        """Test that ontology_term_synonyms returns the correct synonyms when term exists."""
        mock_ontology_parser.get_term_synonyms.return_value = ["synonym1", "synonym2"]

        result = ontology_term_synonyms("CL:0000003")

        self.assertEqual(result, ["synonym1", "synonym2"])
        mock_ontology_parser.get_term_synonyms.assert_called_once_with("CL:0000003")


if __name__ == "__main__":
    unittest.main()
