import json
from unittest.mock import MagicMock, patch

from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


def mock_put_object(bucket_name, key_name, file_content_str):
    assert bucket_name == "cellguide-data-public-dev"
    return None


class TestPostCellGuide(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.test_cellguide_description_upload = json.dumps(
            dict(
                cell_onthology_id="CL_0000030",
                description="this is a description",
                references=["https://doi.org/10.1073/pnas.97.12.6242", "https://f1000research.com/articles/9-233/v1"],
            )
        )

    def test__upload_description__no_auth(self):
        response = self.app.post("/cellguide/v1/upload", self.test_cellguide_description_upload)
        self.assertEqual(401, response.status_code)

    @patch("backend.layers.thirdparty.s3_provider_mock.MockS3Provider")
    def test__upload_description__OK(self, mock_s3_client):
        mock_s3 = MagicMock()
        mock_s3.put_object = mock_put_object
        mock_s3_client.return_value = mock_s3
        headers = self.make_cxg_admin_header()
        response = self.app.post(
            "/cellguide/v1/upload",
            headers=headers,
            data=self.test_cellguide_description_upload,
        )
        self.assertIn("cell_onthology_id", response.json.keys())
        self.assertIn("description", response.json.keys())
        self.assertIn("references", response.json.keys())
        self.assertEqual(201, response.status_code)
