from unittest.mock import patch

from tests.test_utils.mocks import mock_get_marker_gene_data
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


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
        response = self.app.get("/cellguide/v1/marker_genes")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data, self.mock_marker_gene_data)
