import json
import unittest

from backend.cellguide.pipeline.canonical_marker_genes.canonical_markers import CanonicalMarkerGenesCompiler
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


def compare_dicts(dict1, dict2):
    """
    This function recursive compares two dictionaries handling cases where the values
    are unordered arrays with elements that could be dictionaries.
    """
    if len(dict1) != len(dict2):
        return False

    for key in dict1:
        if key not in dict2:
            return False

        value1 = dict1[key]
        value2 = dict2[key]

        if isinstance(value1, dict) and isinstance(value2, dict):
            if not compare_dicts(value1, value2):
                return False
        elif isinstance(value1, list) and isinstance(value2, list):
            if len(value1) != len(value2):
                return False
            # check if the lists contain dictionaries as elements
            if len(value1) > 0 and isinstance(value1[0], dict) and isinstance(value2[0], dict):
                for i in range(len(value1)):
                    if not compare_dicts(value1[i], value2[i]):
                        return False
            elif sorted(value1) != sorted(value2):
                return False
        else:
            if value1 != value2:
                return False

    return True


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
            marker_gene_compiler = CanonicalMarkerGenesCompiler(wmg_tissues, wmg_human_genes)

            canonical_marker_genes = convert_dataclass_to_dict_and_strip_nones(
                marker_gene_compiler.get_processed_asctb_table_entries()
            )

            self.assertTrue(compare_dicts(canonical_marker_genes, expected__canonical_marker_genes))
