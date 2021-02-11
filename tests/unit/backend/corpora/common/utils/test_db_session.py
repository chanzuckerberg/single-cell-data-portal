from backend.corpora.common.corpora_orm import DbDataset
from backend.corpora.common.utils.db_session import db_session_manager, clone
from backend.corpora.common.utils.exceptions import CorporaException
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDBSessionManager(DataPortalTestCase):
    def test_positive(self):
        with self.assertRaises(CorporaException):
            with db_session_manager() as session:
                session.query(DbDataset).filter([DbDataset.id == "test_dataset_id"])[0].update({DbDataset.id: None})


class TestClone(DataPortalTestCase):
    def test_positive(self):
        revision = 1
        dataset_id = "1234"
        dataset = DbDataset(id=dataset_id, revision=revision)
        new_pk = "new_primary_key"
        clone_dataset = clone(dataset, primary_key={"id": new_pk})
        self.assertNotEqual(dataset.id, clone_dataset.id)
        self.assertIsInstance(clone_dataset, dataset.__class__)
        self.assertEqual(clone_dataset.id, new_pk)
        self.assertEqual(clone_dataset.revision, revision)

    def test_bad_field(self):
        revision = 1
        dataset_id = "1234"
        dataset = DbDataset(id=dataset_id, revision=revision)
        new_pk = "new_primary_key"
        with self.assertRaises(TypeError):
            clone(dataset, primary_key={"id": new_pk}, fake_field=False)
