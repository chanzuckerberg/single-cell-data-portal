import unittest

from backend.corpora.common.corpora_orm import DbDataset
from backend.corpora.common.utils.db_utils import DbUtils, db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException


class TestDbUtils(unittest.TestCase):
    def test__single_session(self):
        """Test that only a single sessions is created"""
        db1 = DbUtils()
        db2 = DbUtils()

        self.assertEqual(db1.session, db2.session)


class TestDBSessionManager(unittest.TestCase):
    def test_positive(self):
        with self.assertRaises(CorporaException):
            with db_session_manager() as manager:
                manager.session.query(DbDataset).filter([DbDataset.id == "test_dataset_id"]).update(
                    {DbDataset.id: None}
                )
