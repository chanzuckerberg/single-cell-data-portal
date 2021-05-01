import unittest

from backend.corpora.common.utils.db_session import DBSessionMaker
from tests.unit.backend.fixtures.generate_data_mixin import GenerateDataMixin
from tests.unit.backend.fixtures.test_db import TestDatabaseManager


class DataPortalTestCase(GenerateDataMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        TestDatabaseManager.initialize_db()

    def setUp(self):
        super().setUp()
        self.session = DBSessionMaker().session()

    def tearDown(self) -> None:
        self.session.close()
        super().tearDown()

    @staticmethod
    def reinitialize_database():
        TestDatabaseManager.is_initialized = False
