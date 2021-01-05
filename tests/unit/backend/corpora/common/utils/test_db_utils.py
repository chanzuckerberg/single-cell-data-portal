from backend.corpora.common.corpora_orm import DbDataset
from backend.corpora.common.utils.database_session import DatabaseSession, db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDbUtils(DataPortalTestCase):
    def test__single_session(self):
        """Test that only a single sessions is created"""
        db1 = DatabaseSession()
        db2 = DatabaseSession()

        self.assertEqual(db1.session, db2.session)


class TestDBSessionManager(DataPortalTestCase):
    def test_positive(self):
        with self.assertRaises(CorporaException):
            with db_session_manager() as manager:
                manager.session.query(DbDataset).filter([DbDataset.id == "test_dataset_id"]).update(
                    {DbDataset.id: None}
                )
