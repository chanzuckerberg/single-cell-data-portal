import unittest

from backend.domain.ontology_labels import ontology_term_label, gene_term_label, enrich_development_stage_with_ancestors


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
            "MONDO:0000001": "disease or disorder",
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
            "ENST00000456328": "DDX11L1-202",
            "ENSMUSG00000102693": "4933401J01Rik",
            "ENSMUST00000193812": "4933401J01Rik-201",
            "ENSSASG00005000002": "ORF1ab",
            "ENSSAST00005000002": "ORF1ab",
        }
        for gene_id, expected_gene_label in test_cases.items():
            with self.subTest(gene_id):
                self.assertEqual(gene_term_label(gene_id), expected_gene_label)

    def test__enrich_development_stage_with_ancestors__no_items__expands_correctly(self):
        self.assertEqual(
                [],
                enrich_development_stage_with_ancestors([""]))

    def test__enrich_development_stage_witha_ancestors__single_item__expands_correctly(self):
        self.assertEqual(
                ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
                enrich_development_stage_with_ancestors(["HsapDv:0000008"]))

    def test__enrich_development_stage_with_ancestors__multiple_items__expand_correctly(self):
        self.assertEqual(
                ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001",
                 "HsapDv:0000033",
                 "HsapDv:0000009"],
                enrich_development_stage_with_ancestors(["HsapDv:0000008", "HsapDv:0000033"]))

    def test__enrich_development_stage_with_ancestors__unknown_items__includes_unknown_items(self):
        self.assertEqual(
                ['unknown_term'],
                enrich_development_stage_with_ancestors(["unknown_term"]))


if __name__ == "__main__":
    unittest.main()
