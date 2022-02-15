import json
import unittest
import uuid

from backend.corpora.api_server.app import app
from unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup


class WmgApiTests(unittest.TestCase):

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    def test__latest_snapshot__returns_200(self):
        response = self.app.get("/wmg/v1/snapshots/latest")

        self.assertEqual(200, response.status_code)

    def test__latest_snapshot__returns_valid_response_body(self):
        response = self.app.get("/wmg/v1/snapshots/latest")

        snapshot_id = json.loads(response.data)['snapshot_id']
        try:
            uuid.UUID(snapshot_id)
        except ValueError:
            self.fail()

    def test__primary_filter_dimensions__returns_200(self):
        snapshot_id = uuid.uuid4()

        response = self.app.get(f"/wmg/v1/snapshots/{snapshot_id}/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    def test__primary_filter_dimensions__returns_valid_response_body(self):
        snapshot_id = uuid.uuid4()

        response = self.app.get(f"/wmg/v1/snapshots/{snapshot_id}/primary_filter_dimensions")

        expected = dict(organism_terms=[dict(oid1='olbl1'),
                                        dict(oid2='olbl2')],
                        tissue_type_terms=[dict(ttid1='ttlbl1'),
                                           dict(ttid2='ttlbl2')])

        self.assertEqual(expected, json.loads(response.data))

    def test__query__minimal_valid_request_returns_200(self):
        snapshot_id = uuid.uuid4()

        request = dict(
                filter=dict(
                        gene_term_ids=['gene1'],
                        organism_term_id='organism1',
                        tissue_type_term_ids=['tissuetype1']),
                response_option='include_filter_dims_include_dataset_links')

        response = self.app.post(f"/wmg/v1/snapshots/{snapshot_id}/query", json=request)

        self.assertEqual(200, response.status_code)

    def test__query__minimal_valid_request_returns_valid_response_body(self):
        snapshot_id = uuid.uuid4()

        request = dict(
                filter=dict(
                        gene_term_ids=['gene1'],
                        organism_term_id='organism1',
                        tissue_type_term_ids=['tissuetype1']),
                response_option='include_filter_dims_include_dataset_links')

        response = self.app.post(f"/wmg/v1/snapshots/{snapshot_id}/query", json=request)

        self.assertEqual(dict(expression_summary={},
                              term_id_labels=dict(
                                      genes=[],
                                      cell_types=[])),
                         json.loads(response.data))

    def test__query__invalid_request_returns_400(self):
        snapshot_id = uuid.uuid4()

        response = self.app.post(f"/wmg/v1/snapshots/{snapshot_id}/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property",
                         json.loads(response.data)['detail'])
