import json
import unittest

from backend.cellguide.pipeline.canonical_marker_genes.canonical_markers import CanonicalMarkerGenesCompiler
from backend.cellguide.pipeline.canonical_marker_genes.utils import clean_doi
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class CanonicalMarkerGeneCompilerTests(unittest.TestCase):
    def test__ontology_tree_builder(self):
        with open("tests/unit/cellguide_pipeline/fixtures/canonical_marker_genes.json", "r") as f:
            expected__canonical_marker_genes = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            wmg_tissues = [
                next(iter(i.keys())) for i in snapshot.primary_filter_dimensions["tissue_terms"]["NCBITaxon:9606"]
            ]
            wmg_human_genes = [
                next(iter(i.values())) for i in snapshot.primary_filter_dimensions["gene_terms"]["NCBITaxon:9606"]
            ]
            marker_gene_compiler = CanonicalMarkerGenesCompiler(
                wmg_tissues=wmg_tissues, wmg_human_genes=wmg_human_genes
            )

            canonical_marker_genes = convert_dataclass_to_dict_and_strip_nones(
                marker_gene_compiler.get_processed_asctb_table_entries()
            )

            self.assertTrue(compare_dicts(canonical_marker_genes, expected__canonical_marker_genes))


class CanonicalMarkerGeneCompilerUtilsTests(unittest.TestCase):
    def test__clean_doi(self):

        test_cases = [
            ("10.1016/j.cell.2019.11.025", "10.1016/j.cell.2019.11.025"),
            ("DOI: 10.1016/j.cell.2019.11.025.", "10.1016/j.cell.2019.11.025"),
            (" DOI: 10.1016/j.cell.2019.11.025 ", "10.1016/j.cell.2019.11.025"),
            ("10.1016/j.cell.2019.11.025. ", "10.1016/j.cell.2019.11.025"),
            ("", ""),
            ("DOI: ", ""),
            (" . ", "."),
        ]

        for doi, expected in test_cases:
            with self.subTest(doi=doi):
                self.assertEqual(clean_doi(doi), expected)
