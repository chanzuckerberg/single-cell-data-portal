import unittest

from backend.corpora.common.utils.db_session import DBSessionMaker
from tests.unit.backend.fixtures.generate_data_mixin import GenerateDataMixin
from tests.unit.backend.fixtures.test_db import TestDatabaseManager


class DataPortalTestCase(GenerateDataMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        TestDatabaseManager.initialize_db()

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        self.session = DBSessionMaker().session()

    def tearDown(self) -> None:
        self.session.close()
