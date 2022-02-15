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


if __name__ == '__main__':
    unittest.main()
