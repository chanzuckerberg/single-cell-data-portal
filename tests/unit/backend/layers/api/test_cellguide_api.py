import json
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


class TestPostCellGuide(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.test_cellguide_desciprtion_upload = dict(
            cell_onthology_id="CL_0000030", description="this is a description", contact_name="john doe", references= ["https://doi.org/10.1073/pnas.97.12.6242", "https://f1000research.com/articles/9-233/v1"]
        )

    def test__upload_description__no_auth(self):
        response = self.app.post("/cellguide/v1/upload", data=json.dumps(self.test_cellguide_desciprtion_upload))
        self.assertEqual(401, response.status_code)

    def test__upload_description__OK(self):
        response = self.app.post(
            "/cellguide/v1/upload", headers=self.make_cxg_admin_header(), data=json.dumps(self.test_cellguide_desciprtion_upload)
        )
        self.assertIn("cell_onthology_id", response.json.keys())
        self.assertIn("description", response.json.keys())
        self.assertIn("references", response.json.keys())
        self.assertEqual(201, response.status_code)
