import unittest
from unittest.mock import patch

from backend.api_server.app import app
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class DeAPIV1Tests(unittest.TestCase):
    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    @patch("backend.de.api.v1.load_snapshot")
    def test__differentialExpression_returns_expected_results(self, load_snapshot):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot
            pass
