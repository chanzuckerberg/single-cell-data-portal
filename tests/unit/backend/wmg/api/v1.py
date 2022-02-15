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

    def test__latest_snapshot__returns_correct_body(self):
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

    def test__primary_filter_dimensions__returns_correct_body(self):
        snapshot_id = uuid.uuid4()

        response = self.app.get(f"/wmg/v1/snapshots/{snapshot_id}/primary_filter_dimensions")

        expected = dict(organism_terms=[dict(id='oid1', label='olbl1'),
                                        dict(id='oid2', label='olbl2')],
                        tissue_type_terms=[dict(id='ttid1', label='ttlbl1'),
                                           dict(id='ttid2', label='ttlbl2')])
        self.assertEqual(expected, json.loads(response.data))
