from unittest.mock import patch

from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


def mock_get_marker_gene_data():
    return {
        "Homo sapiens": {
            "brain": {
                "CL:0000540": [
                    {"marker_score": 0.95, "me": 0.5, "pc": 0.1, "gene": "Gene1"},
                    {"marker_score": 0.90, "me": 0.4, "pc": 0.2, "gene": "Gene2"},
                ]
            }
        }
    }


class TestMarkerGenesAPI(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.mock_marker_gene_data = mock_get_marker_gene_data()

    @patch(
        "backend.cellguide.api.v1.marker_genes.actions.get_marker_gene_data",
        new=mock_get_marker_gene_data,
    )
    def test_get_marker_genes(self):
        response = self.app.get(
            "/cellguide/v1/marker_genes?organism=NCBITaxon:9606&tissue=UBERON:0000955&cell_type=CL:0000540"
        )
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0]["gene"], "Gene1")
        self.assertEqual(data[1]["gene"], "Gene2")
