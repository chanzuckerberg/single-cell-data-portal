import json

from moto import mock_aws

from backend.layers.thirdparty.s3_provider import S3Provider
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


class TestPostCellGuide(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.payload = dict(
            cell_ontology_term_id="CL_0000030",
            description="this is a description",
            references=["https://doi.org/10.1073/pnas.97.12.6242", "https://f1000research.com/articles/9-233/v1"],
        )

    def test__upload_description__no_auth(self):
        response = self.app.post(
            "/cellguide/v1/upload",
            headers=self.make_not_auth_header(),
            data=self.payload,
        )
        data = json.loads(response.data.decode())
        self.assertEqual("Unauthorized", data["title"])
        self.assertEqual(401, response.status_code)

    def test__upload_description__wrong_auth(self):
        response = self.app.post(
            "/cellguide/v1/upload",
            headers=self.make_super_curator_header(),
            data=self.payload,
        )
        data = json.loads(response.data.decode())
        self.assertEqual("Bad Request", data["title"])
        self.assertEqual(400, response.status_code)

    @mock_aws
    def test__upload_description__OK(self):
        s3 = S3Provider()
        s3.create_bucket(bucket_name="cellguide-data-public-dev", location="us-west-2")

        response = self.app.post(
            "/cellguide/v1/upload",
            headers=self.make_cxg_admin_header(),
            data=json.dumps(self.payload),
        )

        self.assertIn("cell_ontology_term_id", response.json.keys())
        self.assertIn("description", response.json.keys())
        self.assertIn("references", response.json.keys())
        self.assertEqual(201, response.status_code)

    def test__upload_bad_data(self):
        self.payload["cell_ontology_term_id"] = "CL_0000xx03"
        response = self.app.post(
            "/cellguide/v1/upload",
            headers=self.make_cxg_admin_header(),
            data=json.dumps(self.payload),
        )
        self.assertEqual(403, response.status_code)

        self.payload["cell_ontology_term_id"] = "CL_0000030"
        self.payload["references"] = ["https://notvalidurl"]
        response = self.app.post(
            "/cellguide/v1/upload",
            headers=self.make_cxg_admin_header(),
            data=json.dumps(self.payload),
        )
        self.assertEqual(403, response.status_code)
