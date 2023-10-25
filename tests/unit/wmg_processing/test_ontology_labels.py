import unittest

from backend.wmg.data.ontology_labels import ethnicity_term_label, gene_term_label, ontology_term_label


class OntologyLabelTests(unittest.TestCase):
    def test__ontology_term_label(self):
        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        # This test is fragile, in that changes to our ontologies may break these these tests in the future, but that
        # should happen infrequently
        test_cases = {
            "EFO:0000001": "experimental factor",
            "CL:0000000": "cell",
            "HsapDv:0000000": "human life cycle stage",
            "PATO:0000001": "quality",
            "NCBITaxon:10001": "Sciurus",
            "MmusDv:0000000": "mouse life cycle stage",
            "MONDO:0000001": "disease",
            "UBERON:0002539": "pharyngeal arch",
            "UBERON:0000995 (organoid)": "uterus (organoid)",
            "CL:0000082 (cell culture)": "epithelial cell of lung (cell culture)",
        }
        for ontology_term_id, expected_ontology_term_label in test_cases.items():
            with self.subTest(ontology_term_id):
                self.assertEqual(ontology_term_label(ontology_term_id), expected_ontology_term_label)

    def test__ethnicity_term_label(self):
        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        # This test is fragile, in that changes to our ontologies may break these these tests in the future, but that
        # should happen infrequently
        test_cases = {
            "HANCESTRO:0003": "country",
            "HANCESTRO:0005,HANCESTRO:0014": "European,Hispanic or Latin American",
        }
        for ontology_term_id, expected_ontology_term_label in test_cases.items():
            with self.subTest(ontology_term_id):
                self.assertEqual(ethnicity_term_label(ontology_term_id), expected_ontology_term_label)

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


if __name__ == "__main__":
    unittest.main()
