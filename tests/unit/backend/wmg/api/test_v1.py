import json
import unittest
import uuid
from pprint import pprint

from backend.corpora.api_server.app import app
from backend.wmg.api import v1
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup


class WmgApiV1Tests(unittest.TestCase):
    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    def test__primary_filter_dimensions__returns_200(self):
        response = self.app.get(f"/wmg/v1/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    def test__primary_filter_dimensions__returns_valid_response_body(self):
        response = self.app.get(f"/wmg/v1/primary_filter_dimensions")

        expected = dict(
                snapshot_id=v1.DUMMY_SNAPSHOT_UUID,
                organism_terms=[dict(oid1="olbl1"), dict(oid2="olbl2")],
                tissue_type_terms=[dict(ttid1="ttlbl1"), dict(ttid2="ttlbl2")],
        )

        self.assertEqual(expected, json.loads(response.data))

    def test__query__minimal_valid_request_returns_200(self):
        request = dict(
            filter=dict(gene_term_ids=["gene1"], organism_term_id="organism1", tissue_type_term_ids=["tissuetype1"]),
            response_option="include_filter_dims_include_dataset_links",
        )

        response = self.app.post(f"/wmg/v1/query", json=request)

        self.assertEqual(200, response.status_code)

    def test__query__valid_request_returns_valid_response_body(self):
        request = dict(
            filter=dict(
                gene_term_ids=["gene1", "gene2"],
                organism_term_id="organism1",
                tissue_type_term_ids=["tissuetype1", "tissuetype2"],
            ),
            response_option="include_filter_dims_include_dataset_links",
        )

        response = self.app.post(f"/wmg/v1/query", json=request)

        expected = {
            "snapshot_id": v1.DUMMY_SNAPSHOT_UUID,
            "expression_summary": {
                "gene1": {
                    "tissuetype1": [
                        {"id": "CL00000", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00001", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00002", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                    ],
                    "tissuetype2": [
                        {"id": "CL00000", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00001", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00002", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                    ],
                },
                "gene2": {
                    "tissuetype1": [
                        {"id": "CL00000", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00001", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00002", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                    ],
                    "tissuetype2": [
                        {"id": "CL00000", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00001", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                        {"id": "CL00002", "me": 0.0, "n": 0, "pc": 0.0, "tpc": 0.0},
                    ],
                },
            },
            "term_id_labels": {
                "cell_types": [
                    {"CL00000": "CL00000_label"},
                    {"CL00001": "CL00001_label"},
                    {"CL00002": "CL00002_label"},
                ],
                "genes": [{"gene1": "gene1_label"}, {"gene2": "gene2_label"}],
            },
        }
        pprint(expected)
        self.assertEqual(expected, json.loads(response.data))

    def test__query__invalid_request_returns_400(self):
        response = self.app.post(f"/wmg/v1/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property", json.loads(response.data)["detail"])
