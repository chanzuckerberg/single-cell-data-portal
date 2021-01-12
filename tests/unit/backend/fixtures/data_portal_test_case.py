import unittest

from tests.unit.backend.fixtures.test_db import TestDatabaseManager


class DataPortalTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        TestDatabaseManager.initialize_db()

    @classmethod
    def tearDownClass(cls):
        pass
