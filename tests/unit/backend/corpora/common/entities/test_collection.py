import unittest
from datetime import datetime

from sqlalchemy.exc import SQLAlchemyError

from backend.corpora.common.corpora_orm import CollectionLinkType, DbCollectionLink, CollectionVisibility, DbDataset
from backend.corpora.common.entities import Dataset
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.utils.db_utils import DbUtils
from tests.unit.backend.utils import BogusCollectionParams, BogusDatasetParams


class TestCollection(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_collection_id"
        self.visibility = CollectionVisibility.PUBLIC.name
        self.db = DbUtils()

    def tearDown(self):
        self.db.session.rollback()
        self.db.close()

    def test__get__ok(self):
        key = (self.uuid, self.visibility)

        collection = Collection.get(key)

        # Verify Columns
        self.assertEqual(collection.name, "test_collection")
        self.assertEqual(collection.owner, "test_user_id")

        # Verify Dataset relationship
        dataset = collection.datasets[0]
        self.assertIsInstance(dataset, DbDataset)
        self.assertEqual(dataset.id, "test_dataset_id")
        self.assertEqual(len(dataset.assay), 1)
        self.assertDictEqual(dataset.assay[0], {"ontology_term_id": "test_obo", "label": "test_assay"})

        # Verify Link relationship
        self.assertIsInstance(collection.links[0], DbCollectionLink)
        self.assertEqual(collection.links[0].id, "test_collection_link_id")

    def test__get__does_not_exist(self):
        non_existent_key = ("non_existent_id", self.visibility)

        self.assertEqual(Collection.get(non_existent_key), None)

    def test__get__invalid_visibility(self):
        invalid_visibility_key = (self.uuid, "invalid_visibility")
        with self.assertRaises(SQLAlchemyError):
            Collection.get(invalid_visibility_key)

    def test__create__ok(self):
        """
        Create a collection with a variable number of links.
        """

        link_params = {"link_url": "fake_url", "link_type": CollectionLinkType.PROTOCOL.name}
        collection_params = BogusCollectionParams.get()

        for i in range(3):
            with self.subTest(i):
                collection = Collection.create(links=[link_params] * i, **collection_params)

                collection_key = (collection.id, collection.visibility)
                expected_links = collection.links

                # Expire all local object and retrieve them from the DB to make sure the transactions went through.
                Collection.db.session.expire_all()

                actual_collection = Collection.get(collection_key)
                self.assertEqual(collection_key, (actual_collection.id, actual_collection.visibility))
                self.assertCountEqual(expected_links, actual_collection.links)

    def test__list_in_time_range__ok(self):
        created_before = Collection.create(**BogusCollectionParams.get(), created_at=datetime.fromtimestamp(10))
        from_date = 20
        created_inbetween = Collection.create(**BogusCollectionParams.get(), created_at=datetime.fromtimestamp(30))
        to_date = 40
        created_after = Collection.create(**BogusCollectionParams.get(), created_at=datetime.fromtimestamp(50))

        with self.subTest("from_date"):
            # Collections from_date are returned.
            actual_collections = Collection.list_attributes_in_time_range(from_date=from_date)
            self.assertTrue(all([p["created_at"].timestamp() > from_date for p in actual_collections]))
            expected_ids = [created_inbetween.id, created_after.id, "test_collection_id"]
            actual_ids = [p["id"] for p in actual_collections]
            # Check if the test ids we created are present.
            # As a result of other tests, more collections have likely been created and will be return in the results,
            # so we can't do an exact match.
            self.assertTrue(set(expected_ids).issubset(actual_ids))

        with self.subTest("to_date"):
            # Collections to_date are returned.
            actual_collections = Collection.list_attributes_in_time_range(to_date=to_date)
            self.assertTrue(all([p["created_at"].timestamp() < to_date for p in actual_collections]))
            expected_ids = [created_before.id, created_inbetween.id]
            actual_ids = [p["id"] for p in actual_collections]
            self.assertCountEqual(expected_ids, actual_ids)

        with self.subTest("from_date->to_date"):
            # Collections between to_date and from_date are returned.
            actual_collections = Collection.list_attributes_in_time_range(to_date=to_date, from_date=from_date)
            self.assertTrue(all([p["created_at"].timestamp() > from_date for p in actual_collections]))
            self.assertTrue(all([p["created_at"].timestamp() < to_date for p in actual_collections]))
            expected_ids = [created_inbetween.id]
            actual_ids = [p["id"] for p in actual_collections]
            self.assertCountEqual(expected_ids, actual_ids)

        with self.subTest("No parameters"):
            """All collections are returned."""
            actual_collections = Collection.list_attributes_in_time_range()
            expected_ids = [created_before.id, created_inbetween.id, created_after.id]
            actual_ids = [p["id"] for p in actual_collections]
            # Check if the test ids we created are present.
            # As a result of other tests, more collections have likely been created and will be return in the results,
            # so we can't do an exact match.
            self.assertTrue(set(expected_ids).issubset(actual_ids))

    def test__list_collections_in_time_range___ok(self):
        """Public collections are returned"""
        from_date = 10
        created_at = datetime.fromtimestamp(20)
        to_date = 30
        expected_collection = Collection.create(
            **BogusCollectionParams.get(created_at=created_at, visibility=CollectionVisibility.PUBLIC.name)
        )

        # Generate a collection with Private visibility. This should not show up in the results.
        Collection.create(
            **BogusCollectionParams.get(created_at=created_at, visibility=CollectionVisibility.PRIVATE.name)
        )

        actual_collections = Collection.list_collections_in_time_range(to_date=to_date, from_date=from_date)
        expected_collections = [{"created_at": created_at, "id": expected_collection.id}]
        self.assertCountEqual(expected_collections, actual_collections)

    def test__list__ok(self):
        generate = 2
        generated_ids = [Collection.create(**BogusCollectionParams.get()).id for _ in range(generate)]
        collections = Collection.list()
        self.assertTrue(set(generated_ids).issubset([p.id for p in collections]))

    def test__cascade_delete_collection__ok(self):
        test_collection = Collection.create(**BogusCollectionParams.get(links=[{}]))
        db = DbUtils()
        test_collection_id = test_collection.id
        test_link_id = test_collection.links[0].id

        # Check if the link exists.
        expected_id = test_link_id
        actual_id = db.query([DbCollectionLink], [DbCollectionLink.id == test_link_id])[0].id
        self.assertEqual(expected_id, actual_id)

        # The Collection is deleted
        db.session.delete(test_collection.db_object)
        expected_results = None
        actual_results = Collection.get_collection(test_collection_id)
        self.assertEqual(expected_results, actual_results)

        # The link should also be deleted.
        expected_results = []
        actual_results = db.query([DbCollectionLink], [DbCollectionLink.id == test_link_id])
        self.assertEqual(expected_results, actual_results)

    def test__cascade_delete_collection_with_dataset__ok(self):
        db = DbUtils()
        test_collection = Collection.create(**BogusCollectionParams.get())
        expected_collection_id = test_collection.id
        test_dataset = Dataset.create(
            **BogusDatasetParams.get(collection_id=test_collection.id, collection_visibility=test_collection.visibility)
        )
        expected_dataset_id = test_dataset.id

        # The Collection is deleted
        db.session.delete(test_collection.db_object)
        db.session.expire_all()
        expected_results = None
        actual_results = Collection.get_collection(expected_collection_id)
        self.assertEqual(expected_results, actual_results)

        # The dataset is delete
        expected_results = None
        actual_results = Dataset.get(expected_dataset_id)
        self.assertEqual(expected_results, actual_results)
