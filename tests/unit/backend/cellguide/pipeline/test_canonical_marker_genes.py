import json
import unittest
from unittest.mock import Mock, patch

from backend.cellguide.pipeline.canonical_marker_genes.canonical_markers import CanonicalMarkerGenesCompiler
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils import compare_dicts
from tests.test_utils.mocks import mock_get_asctb_master_sheet, mock_get_title_and_citation_from_doi
from tests.unit.backend.cellguide.pipeline.constants import (
    CANONICAL_MARKER_GENES_FIXTURE_FILENAME,
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class CanonicalMarkerGeneCompilerTests(unittest.TestCase):
    def test__canonical_marker_genes(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CANONICAL_MARKER_GENES_FIXTURE_FILENAME}", "r") as f:
            expected__canonical_marker_genes = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            wmg_tissues = [
                next(iter(i.keys())) for i in snapshot.primary_filter_dimensions["tissue_terms"]["NCBITaxon:9606"]
            ]
            wmg_human_genes = [
                next(iter(i.values())) for i in snapshot.primary_filter_dimensions["gene_terms"]["NCBITaxon:9606"]
            ]

            with (
                patch(
                    "backend.cellguide.pipeline.canonical_marker_genes.canonical_markers.get_asctb_master_sheet",
                    new=mock_get_asctb_master_sheet,
                ),
                patch(
                    "backend.cellguide.pipeline.canonical_marker_genes.canonical_markers.CrossrefProvider",
                    new=Mock(get_title_and_citation_from_doi=mock_get_title_and_citation_from_doi),
                ),
            ):
                marker_gene_compiler = CanonicalMarkerGenesCompiler(
                    wmg_tissues=wmg_tissues, wmg_human_genes=wmg_human_genes
                )

                canonical_marker_genes = convert_dataclass_to_dict_and_strip_nones(
                    marker_gene_compiler.get_processed_asctb_table_entries()
                )
            self.assertTrue(compare_dicts(canonical_marker_genes, expected__canonical_marker_genes))
