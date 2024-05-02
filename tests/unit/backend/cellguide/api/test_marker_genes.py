from unittest.mock import patch

from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest
from tests.test_utils.mocks import mock_get_marker_gene_data


class TestMarkerGenesAPI(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.mock_marker_gene_data = mock_get_marker_gene_data()

    @patch(
        "backend.cellguide.api.v1.marker_genes.actions.get_marker_gene_data",
        new=mock_get_marker_gene_data,
    )
    def test_get_marker_genes(self):
        # Test case: Organism, tissue, and cell type provided
        response = self.app.get(
            "/cellguide/v1/marker_genes?organism=NCBITaxon:9606&tissue=UBERON:0000955&cell_type=CL:0000540"
        )
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data, self.mock_marker_gene_data["Homo sapiens"]["brain"]["CL:0000540"])

        # Test case: Only organism is provided
        response = self.app.get("/cellguide/v1/marker_genes?organism=NCBITaxon:9606")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        expected_data = self.mock_marker_gene_data["Homo sapiens"]
        self.assertEqual(data, expected_data)

        # Test case: Organism and tissue is provided
        response = self.app.get("/cellguide/v1/marker_genes?organism=NCBITaxon:9606&tissue=brain")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        expected_data = self.mock_marker_gene_data["Homo sapiens"]["brain"]
        self.assertEqual(data, expected_data)

        # Test case: Organism and cell type is provided, agnostic to tissues
        response = self.app.get("/cellguide/v1/marker_genes?organism=NCBITaxon:9606&cell_type=CL:0000540")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        expected_data = self.mock_marker_gene_data["Homo sapiens"]["All Tissues"]["CL:0000540"]
        self.assertEqual(data, expected_data)
