import unittest

from tests.unit.backend.fixtures.test_db import TestDatabaseManager


class DataPortalTestCase(unittest.TestCase):
    def setUp(self):
        TestDatabaseManager.initialize_db()

    def tearDown(self):
        pass
