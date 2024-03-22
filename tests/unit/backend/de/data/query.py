import unittest

from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"


class DeQueryTest(unittest.TestCase):
    def test__correct_cube_selected_based_on_query(self):
        with load_realistic_test_snapshot(TEST_SNAPSHOT):
            pass
