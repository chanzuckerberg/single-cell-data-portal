from datetime import datetime
import unittest
import copy
from unittest.mock import Mock
from backend.corpora.common.providers.crossref_provider import CrossrefException, CrossrefProvider
from backend.layers.business.business import BusinessLogic, CollectionQueryFilter, UserInfo
from backend.layers.common.entities import CollectionMetadata, CollectionVersion, DatasetMetadata, DatasetStatus, Link
from backend.layers.persistence.persistence import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock

class BaseBusinessLogicTestCase(unittest.TestCase):

    sample_collection_metadata = CollectionMetadata(
        "test collection 1", 
        "description of test collection 1",
        "test_user_1",
        "scientist",
        "scientist@czi.com",
        []
    )

    sample_dataset_metadata = DatasetMetadata(
        "test_organism",
        "test_tissue",
        "test_assay",
        "test_disease", 
        "test_sex",
        "test_self_reported_ethnicity",
        "test_development_stage",
        "test_cell_type",
        10
    )

    test_user_name = "test_user_1"

    user_info = UserInfo(
        "test_user_1",
        "fake token"
    )

    def setUp(self) -> None:
        self.database_provider = DatabaseProviderInterface() # replace with DatabaseProviderMock()

        # By default does nothing. Can be mocked by single test cases.
        self.crossref_provider = CrossrefProviderInterface()
        self.step_function_provider = StepFunctionProviderInterface()

        self.business_logic = BusinessLogic(
            database_provider=self.database_provider, 
            crossref_provider=self.crossref_provider,
            step_function_provider=self.step_function_provider,
        )
        return super().setUp()

    def initialize_empty_private_collection(self, owner: str = test_user_name, metadata = sample_collection_metadata) -> CollectionVersion:
        """
        Initializes a private collection to be used for testing, with no datasets
        """
        version = self.database_provider.initialize_canonical_collection(
            owner,
            metadata,
        )
        return version

    def initialize_private_collection(self, owner: str = test_user_name, metadata = sample_collection_metadata) -> CollectionVersion:
        """
        Initializes a private collection to be used for testing, with a single dataset
        """
        version = self.initialize_empty_private_collection()
        dataset_version = self.database_provider.initialize_canonical_dataset(
            version.version_id, 
            dataset_metadata=self.sample_dataset_metadata
        )
        self.database_provider.add_dataset_to_collection_version_mapping(
            version.version_id, 
            dataset_version.version_id
        )
        return version
        

    def initialize_public_collection(self, owner: str = test_user_name, published_at: datetime = None, metadata = sample_collection_metadata) -> CollectionVersion:
        """
        Initializes a public collection to be used for testing, with a single dataset
        """
        version = self.initialize_private_collection(owner, metadata)
        self.database_provider.finalize_collection_version(version.collection_id, version.version_id, published_at)
        return version
        


class TestCreateCollection(BaseBusinessLogicTestCase):

    def create_collection_ok(self):
        """
        A collection can be created using `create_collection`
        """
        cv = self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)
        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(cv, cv_from_database)

    def create_collection_with_links_ok(self):
        """
        A collection with links can be created using `create_collection`
        """
        good_links = [
            Link("test link 1", "protocol", "http://example.com/protocol"),
            Link("test link 2", "other", "http://example.com/other"),
        ]
        self.sample_collection_metadata.links = good_links
        cv = self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)
        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(cv.metadata.links, cv_from_database.metadata.links)

    def create_collection_with_bad_links_fail(self):
        """
        A collection where a link isn't a http URI cannot be created
        """
        bad_links = [
            Link("test bad link", "other", "incorrect_url")
        ]
        self.sample_collection_metadata.links = bad_links
        # TODO: Create InvalidLinkException
        # TODO: the current method has error message validation - implement a more sophisticate version
        with self.assertRaises(InvalidLinkException):
            self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)

    def create_collection_with_valid_doi_ok(self):
        """
        A collection can be created with a valid DOI, and the publisher metadata will be added
        to the collection
        """
        links_with_doi = [
            Link("test doi", "doi", "http://good.doi")
        ]
        self.sample_collection_metadata.links = links_with_doi

        expected_publiser_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publiser_metadata)

        cv = self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)

        self.crossref_provider.fetch_metadata.assert_called_with("http://good.doi")

        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(1, len(cv_from_database.metadata.links))
        self.assertEqual(cv_from_database.metadata.links[0].uri, "http://good.doi")
        self.assertIsNotNone(cv_from_database.publisher_metadata)
        self.assertEqual(cv_from_database.publisher_metadata, expected_publiser_metadata)


    def create_collection_with_invalid_doi_fail(self):
        """
        A collection with an invalid DOI will not be created
        """
        links_with_doi = [
            Link("test doi", "doi", "http://bad.doi")
        ]

        self.sample_collection_metadata.links = links_with_doi

        # TODO: make sure that we don't need different actions depending on which exception
        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefException("Error!"))

        cv = self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)
        # TODO: create CollectionCreationException
        with self.assertRaises(CollectionCreationException):
            self.business_logic.create_collection(self.sample_collection_metadata, self.user_info)

    def create_collection_unauthorized_fail(self):
        # TODO: AUTHORIZATION
        pass

class TestDeleteCollection(BaseBusinessLogicTestCase):

    def delete_collection_ok(self):
        # TODO: verify how the current support for collection deletion is
        pass

    def delete_collection_unauthorized_fail(self):
        # TODO: AUTHORIZATION
        pass

class TestGetCollection(BaseBusinessLogicTestCase):

    def get_collection_published_ok(self):
        """
        A published collection can be obtained using `get_collection`
        """
        version = self.initialize_public_collection()
        fetched_version = self.business_logic.get_collection(version.collection_id, self.user_info)
        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

    def get_collection_unpublished_fail(self):
        """
        An unpublished collection should not be obtained using get_collection
        Instead, `get_collection_version` should be used
        """
        version = self.initialize_private_collection()
        fetched_version = self.business_logic.get_collection(version.collection_id, self.user_info)
        self.assertIsNone(fetched_version)

    def get_collection_version_ok(self):
        """
        A collection version can be obtained using `get_collection_version`
        """
        version = self.initialize_private_collection()
        fetched_version = self.business_logic.get_collection_version(version.version_id, self.user_info)
        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

class TestGetAllCollections(BaseBusinessLogicTestCase):

    def setUp(self) -> None:
        self.collection1 = self.initialize_private_collection(owner="test_user_1")
        self.collection2 = self.initialize_private_collection(owner="test_user_2")
        self.collection3 = self.initialize_public_collection(published_at=datetime(2022, 1, 1))
        self.collection4 = self.initialize_public_collection(published_at=datetime(2021, 1, 1))
        self.collection5 = self.initialize_public_collection(published_at=datetime(2020, 1, 1))
        return super().setUp()

    def get_all_collections_ok(self):
        """
        All the collection versions should be returned by `get_collections`.
        """
        # TODO: this method should NOT be used without at least one filter. Maybe add an assertion to block it?
        filter = CollectionQueryFilter(
            from_date = None,
            to_date = None,
            visibility=None,
            owner = None
        )
        versions = self.business_logic.get_collections(filter, self.user_info)
        self.assertEqual(5, len(list(versions)))

    def get_all_collections_with_owner_ok(self):
        """
        Only collection versions with the right owner should be returned
        """
        filter = CollectionQueryFilter(
            from_date = None,
            to_date = None,
            visibility=None,
            owner = "test_user_2"
        )
        # TODO: I believe that `owner` should be specified in user_info. Investigate
        versions = list(self.business_logic.get_collections(filter, self.user_info))
        self.assertEqual(1, len(versions))
        self.assertEqual(versions[0].version_id, self.collection2.version_id)
        

    def get_all_collections_with_date_range_ok(self):
        """
        Only collection versions where `published_at` is within the date range should be returned
        """
        filter = CollectionQueryFilter(
            from_date = "2020-12-01",
            to_date = "2021-12-01",
            visibility=None,
            owner = None
        )
        versions = list(self.business_logic.get_collections(filter, self.user_info))
        self.assertEqual(1, len(versions))
        self.assertEqual(versions[0].version_id, self.collection4.version_id)

    def get_all_collections_visibility_ok(self):
        """
        If visibility = "PRIVATE", only unpublished collections should be returned
        If visibility = "PUBLIC", only published collections should be returned
        """
        filter = CollectionQueryFilter(
            from_date = None,
            to_date = None,
            visibility="PRIVATE",
            owner = None
        )
        versions = list(self.business_logic.get_collections(filter, self.user_info))
        self.assertEqual(2, len(versions))
        for version in versions:
            self.assertIsNone(version.published_at)

        filter = CollectionQueryFilter(
            from_date = None,
            to_date = None,
            visibility="PUBLIC",
            owner = None
        )
        versions = list(self.business_logic.get_collections(filter, self.user_info))
        self.assertEqual(3, len(versions))
        for version in versions:
            self.assertIsNotNone(version.published_at)


class TestUpdateCollection(BaseBusinessLogicTestCase):

    def update_collection_ok(self):
        """
        A collection version should be updated when using `update_collection`
        """
        version = self.initialize_private_collection()
        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
        }

        # TODO: can `owner` be modified by any role?

        self.business_logic.update_collection_version(version.version_id, body, self.user_info)
        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.metadata.name, body["name"])
        self.assertEqual(updated_version.metadata.description, body["description"])
        self.assertEqual(updated_version.metadata.contact_name, body["contact_name"])
        self.assertEqual(updated_version.metadata.contact_email, body["contact_email"])
        self.assertEqual(updated_version.metadata.owner, self.sample_collection_metadata.owner)

    def update_collection_wrong_params_fail(self):
        """
        Updating a collection with a payload with missing fields should fail
        """
        version = self.initialize_private_collection()
        body = {
            "name": "new collection name",
        }
        # TODO: Add more update failure cases?
        # TODO: Check the returned error messages
        with self.assertRaises(CollectionUpdateException):
            self.business_logic.update_collection_version(version.version_id, body, self.user_info)

    def update_published_collection_fail(self):
        """
        Updating a collection version that is published should fail
        """
        version = self.initialize_public_collection()
        # TODO: Create CollectionUpdateException
        with self.assertRaises(CollectionUpdateException):
            self.business_logic.update_collection_version(version.version_id, {"name": "test"}, self.user_info)

    def update_collection_same_doi(self):
        """
        A collection updated with the same DOI should not trigger a Crossref call
        """
        metadata = self.sample_collection_metadata
        links = [Link("test doi", "DOI", "http://test.doi")]
        metadata.links = links

        expected_publiser_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publiser_metadata)

        version = self.initialize_private_collection(metadata=metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
            "links": links
        }

        self.business_logic.update_collection_version(version.version_id, body, self.user_info)
        self.crossref_provider.fetch_metadata.assert_not_called()

        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertIsNotNone(updated_version.publisher_metadata)

        

    def update_collection_change_doi(self):
        """
        A collection updated with a new DOI should get new publisher metadata from Crossref
        """
        metadata = self.sample_collection_metadata
        links = [Link("test doi", "DOI", "http://test.doi")]
        metadata.links = links

        self.crossref_provider.fetch_metadata = Mock(return_value={"authors": ["Test Author"]})

        version = self.initialize_private_collection(metadata=metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
            "links": [Link("new test doi", "DOI", "http://new.test.doi")]
        }

        expected_updated_publisher_metadata = {"authors": ["New Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_updated_publisher_metadata)

        self.business_logic.update_collection_version(version.version_id, body, self.user_info)
        self.crossref_provider.fetch_metadata.assert_called_once()

        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.publisher_metadata, expected_updated_publisher_metadata)


    def update_collection_unauthorized_fail(self):
        """
        An unauthorized attempt to update a collection should be rejected
        """
        # TODO: AUTHORIZATION
        pass

class TestAddDataset(BaseBusinessLogicTestCase):

    def add_dataset_to_collection_ok(self):
        """
        A dataset can be added to a collection when `ingest_dataset` is called.
        The resulting dataset should be empty and in a state ready for processing.
        """
        version = self.initialize_empty_private_collection()
        url = "http://test/dataset.url"
        new_dataset_version_id = self.business_logic.ingest_dataset(version.version_id, url, None, self.user_info)
        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.processing_status, DatasetStatus.WAITING)

    def add_dataset_to_public_collection_fail(self):
        """
        Adding a dataset to a public collection should result in a failure
        """
        version = self.initialize_public_collection()
        url = "http://test/dataset.url"
        with self.assertRaises(DatasetIngestException):
            self.business_logic.ingest_dataset(version.version_id, url, None, self.user_info)

class TestUpdateRemoveDataset(BaseBusinessLogicTestCase):

    def remove_dataset_from_collection_ok(self):
        """
        A dataset can be removed from a collection version using `delete_dataset`
        """
        version = self.initialize_private_collection()
        self.assertEqual(1, len(version.datasets))
        self.business_logic.delete_dataset(version.version_id)

    def remove_dataset_from_public_collection_fail(self):
        """
        Removing a dataset from a published collection should fail
        """
        pass

    def replace_dataset_in_collection_ok(self):
        """
        A dataset can be replaced from a collection version by calling `ingest_dataset`
        and specifying an existing dataset_version_id
        """
        pass

    def replace_dataset_in_public_collection_fail(self):
        """
        Replacing a dataset that belongs to a published collection should fail
        """
        pass

class TestGetDataset(BaseBusinessLogicTestCase):

    def get_all_datasets_ok(self):
        """
        All dataset that belong to a published collection can be retrieved with `get_all_datasets`
        """
        pass

    def get_dataset_artifacts_ok(self):
        """
        Artifacts belonging to a dataset can be obtained with `get_dataset_artifacts`
        """
        pass

class TestDatasetStatus(BaseBusinessLogicTestCase):

    def update_dataset_status_ok(self):
        pass

    def get_dataset_status_ok(self):
        pass

class TestNewCollectionVersion(BaseBusinessLogicTestCase):

    def create_collection_version_ok(self):
        pass

    def create_collection_version_unauthorized_fail(self):
        pass

    def delete_collection_version_ok(self):
        pass

    def delete_collection_version_unauthorized_fail(self):
        pass

    def publish_version_ok(self):
        pass

    def publish_version_with_updated_metadata_ok(self):
        pass

    def publish_version_with_removed_dataset_ok(self):
        pass

    def publish_version_with_added_dataset_ok(self):
        pass

    def publish_version_with_updated_dataset_ok(self):
        pass

    def publish_version_unauthorized_fail(self):
        pass


if __name__ == '__main__':
    unittest.main()

