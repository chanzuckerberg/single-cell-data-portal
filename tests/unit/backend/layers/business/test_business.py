from audioop import add
import unittest
from datetime import datetime
from unittest.mock import Mock

from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface

from backend.corpora.common.providers.crossref_provider import CrossrefException
from backend.layers.business.business import BusinessLogic, CollectionQueryFilter, DatasetArtifactDownloadData
from backend.layers.business.exceptions import CollectionUpdateException, InvalidLinkException, \
    CollectionCreationException, DatasetIngestException, CollectionPublishException
from backend.layers.common.entities import CollectionMetadata, CollectionVersion, DatasetArtifact, DatasetMetadata, DatasetStatus, DatasetVersionId, Link
from backend.layers.persistence.persistence import DatabaseProviderInterface
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock


class BaseBusinessLogicTestCase(unittest.TestCase):

    sample_collection_metadata = CollectionMetadata(
        "test collection 1", 
        "description of test collection 1",
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

    def setUp(self) -> None:
        self.database_provider = DatabaseProviderMock()

        # By default does nothing. Can be mocked by single test cases.
        self.crossref_provider = CrossrefProviderInterface()
        self.step_function_provider = StepFunctionProviderInterface()

        self.business_logic = BusinessLogic(
            database_provider=self.database_provider, 
            crossref_provider=self.crossref_provider,
            step_function_provider=self.step_function_provider,
        )

    def initialize_empty_unpublished_collection(self, owner: str = test_user_name, metadata = sample_collection_metadata) -> CollectionVersion:
        """
        Initializes an unpublished collection to be used for testing, with no datasets
        """
        version = self.database_provider.create_canonical_collection(
            owner,
            metadata,
        )
        return version

    def initialize_unpublished_collection(self, 
        owner: str = test_user_name, 
        metadata = sample_collection_metadata, 
        complete_dataset_ingestion: bool = True
    ) -> CollectionVersion:
        """
        Initializes an unpublished collection to be used for testing, with two datasets.
        By default also completes dataset ingestion (normally, a process that would be done asynchonously).
        Pass `complete_dataset_ingestion=False` if you want to initialize datasets only.
        """
        version = self.initialize_empty_unpublished_collection(owner, metadata)
        for i in range(2):
            dataset_version = self.database_provider.create_canonical_dataset(
                version.version_id, 
                dataset_metadata=self.sample_dataset_metadata
            )
            self.database_provider.add_dataset_to_collection_version_mapping(
                version.version_id, 
                dataset_version.version_id
            )
            if complete_dataset_ingestion:
                self.complete_dataset_processing_with_success(dataset_version.version_id)
        return self.database_provider.get_collection_version(version.version_id)
        
    def initialize_published_collection(self, owner: str = test_user_name, published_at: datetime = datetime.utcnow(), metadata = sample_collection_metadata) -> CollectionVersion:
        """
        Initializes a published collection to be used for testing, with a single dataset
        """
        version = self.initialize_unpublished_collection(owner, metadata)
        self.database_provider.finalize_collection_version(version.collection_id, version.version_id, published_at)
        return self.database_provider.get_collection_version(version.version_id)

    def complete_dataset_processing_with_success(self, dataset_version_id: DatasetVersionId) -> None:
        """
        Test method that "completes" a dataset processing. This is necessary since dataset ingestion
        is a complex process which happens asynchronously, and cannot be easily mocked.
        """
        self.database_provider.add_dataset_artifact(dataset_version_id, "H5AD", "s3://fake-bucket/artifact.h5ad")
        self.database_provider.add_dataset_artifact(dataset_version_id, "CXG", "s3://fake-bucket/artifact.cxg")
        self.database_provider.add_dataset_artifact(dataset_version_id, "RDS", "s3://fake-bucket/artifact.rds")
        self.database_provider.update_dataset_processing_status(dataset_version_id, DatasetStatus.UPLOADED)



class TestCreateCollection(BaseBusinessLogicTestCase):

    def test_create_collection_ok(self):
        """
        A collection can be created using `create_collection`
        """
        collection = self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)
        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(collection, collection_from_database)

    def test_create_collection_with_links_ok(self):
        """
        A collection with links can be created using `create_collection`
        """
        good_links = [
            Link("test link 1", "protocol", "http://example.com/protocol"),
            Link("test link 2", "other", "http://example.com/other"),
        ]
        self.sample_collection_metadata.links = good_links
        collection = self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)
        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(good_links, collection_from_database.metadata.links)

    def test_create_collection_with_bad_links_fail(self):
        """
        A collection where a link isn't a http URI cannot be created
        """
        bad_links = [
            Link("test bad link", "other", "incorrect_url")
        ]
        self.sample_collection_metadata.links = bad_links

        with self.assertRaises(InvalidLinkException) as ex:
            self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)

        self.assertEqual(ex.exception.errors, [
            {"name": f"links[0]", "reason": "Invalid URL.", "value": "incorrect_url"}
        ])

    def test_create_collection_with_valid_doi_ok(self):
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

        collection = self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)

        self.crossref_provider.fetch_metadata.assert_called_with("http://good.doi")

        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(1, len(collection_from_database.metadata.links))
        self.assertEqual(collection_from_database.metadata.links[0].uri, "http://good.doi")
        self.assertIsNotNone(collection_from_database.publisher_metadata)
        self.assertEqual(collection_from_database.publisher_metadata, expected_publiser_metadata)

    def test_create_collection_with_invalid_doi_fail(self):
        """
        A collection with an invalid DOI will not be created
        """
        links_with_doi = [
            Link("test doi", "doi", "http://bad.doi")
        ]

        self.sample_collection_metadata.links = links_with_doi

        # TODO: make sure that we don't need different actions depending on which exception
        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefException("Error!"))

        self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)
        with self.assertRaises(CollectionCreationException):
            self.business_logic.create_collection(self.test_user_name, self.sample_collection_metadata)
        
class TestGetCollectionVersion(BaseBusinessLogicTestCase):

    def test_get_published_collection_version_for_published_collection_ok(self):
        """
        A published collection can be obtained using `get_collection`
        """
        version = self.initialize_published_collection()

        fetched_version = self.business_logic.get_published_collection_version(version.collection_id)

        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

    def test_get_published_collection_version_for_published_collection_is_none(self):
        """
        An unpublished collection should not be obtained using get_collection
        Instead, `get_collection_version` should be used
        """
        version = self.initialize_unpublished_collection()
        fetched_version = self.business_logic.get_published_collection_version(version.collection_id)
        self.assertIsNone(fetched_version)

    def test_get_collection_version_for_unpublished_collection_ok(self):
        """
        An unpublished collection version can be obtained using `get_collection_version`
        """
        version = self.initialize_unpublished_collection()

        fetched_version = self.business_logic.get_collection_version(version.version_id)

        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

    def test_get_collection_version_for_published_collection_ok(self):
        """
        A published collection version can be obtained using `get_collection_version`
        """
        version = self.initialize_published_collection()

        fetched_version = self.business_logic.get_collection_version(version.version_id)

        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

class TestGetAllCollections(BaseBusinessLogicTestCase):

    def test_get_all_collections_unfiltered_ok(self):
        """
        All the collection versions should be returned by `get_collections`, including published, unpublished, and all owners
        """
        self.initialize_unpublished_collection()
        self.initialize_unpublished_collection(owner="test_user_2")
        self.initialize_published_collection()
        self.initialize_published_collection(owner="test_user_2")

        # TODO: this method should NOT be used without at least one filter. Maybe add an assertion to block it?
        filter = CollectionQueryFilter()
        versions = self.business_logic.get_collections(filter)

        self.assertEqual(4, len(list(versions)))

    def test_get_all_collections_with_owner_ok(self):
        """
        Only collection versions with the specified owner should be returned, including both published and unpublished
        """
        self.initialize_unpublished_collection(owner="test_user_1")
        self.initialize_published_collection(owner="test_user_1")
        self.initialize_unpublished_collection(owner="test_user_2")

        filter = CollectionQueryFilter(owner="test_user_1")
        versions = list(self.business_logic.get_collections(filter))

        self.assertEqual(2, len(versions))
        for version in versions:
            self.assertEqual("test_user_1", version.owner)

    def test_get_all_collections_published_ok(self):
        """
        If published filter flag is True, only published collections should be returned
        """
        self.initialize_unpublished_collection()
        self.initialize_published_collection()
        self.initialize_published_collection()

        filter = CollectionQueryFilter(is_published=True)
        versions = list(self.business_logic.get_collections(filter))

        self.assertEqual(2, len(versions))
        for version in versions:
            self.assertIsNotNone(version.published_at)

    def test_get_all_collections_unpublished_ok(self):
        """
        If published filter flag is False, only unpublished collections should be returned
        """
        self.initialize_unpublished_collection()
        self.initialize_unpublished_collection()
        self.initialize_published_collection()

        filter = CollectionQueryFilter(is_published=False)
        versions = list(self.business_logic.get_collections(filter))

        self.assertEqual(2, len(versions))
        for version in versions:
            self.assertIsNone(version.published_at)

        filter = CollectionQueryFilter(is_published=True)

        versions = list(self.business_logic.get_collections(filter))
        self.assertEqual(2, len(versions))
        for version in versions:
            self.assertIsNotNone(version.published_at)


class TestUpdateCollection(BaseBusinessLogicTestCase):
    """
    Tests operations that can update an unpublished collection version. Also tests that these operations cannot be
    performed on a published collection.
    """

    def test_update_collection_ok(self):
        """
        A collection version should be updated when using `update_collection`
        """
        version = self.initialize_unpublished_collection()
        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
        }

        self.business_logic.update_collection_version(version.version_id, body)

        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.metadata.name, body["name"])
        self.assertEqual(updated_version.metadata.description, body["description"])
        self.assertEqual(updated_version.metadata.contact_name, body["contact_name"])
        self.assertEqual(updated_version.metadata.contact_email, body["contact_email"])

    def test_update_collection_wrong_params_fail(self):
        """
        Updating a collection with a payload with missing fields should fail
        """
        version = self.initialize_unpublished_collection()
        body = {
            "name": "new collection name",
        }

        with self.assertRaises(CollectionUpdateException) as ex:
            self.business_logic.update_collection_version(version.version_id, body)

        self.assertEqual(ex.exception.errors, [
            {"name": "description", "reason": "Cannot be blank."},
            {"name": "contact_name", "reason": "Cannot be blank."},
            {"name": "contact_email", "reason": "Cannot be blank."},
        ])

        # TODO: I would recommend adding more validation logic testing (see existing TestVerifyCollection)
        # in a separate unit test case later

    def test_update_published_collection_fail(self):
        """
        Updating a collection version that is published should fail
        """
        version = self.initialize_published_collection()

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.update_collection_version(version.version_id, {"name": "test"})

    def test_update_collection_same_doi(self):
        """
        A collection updated with the same DOI should not trigger a Crossref call
        """
        metadata = self.sample_collection_metadata
        links = [Link("test doi", "DOI", "http://test.doi")]
        metadata.links = links

        expected_publiser_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publiser_metadata)

        version = self.initialize_unpublished_collection(metadata=metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
            "links": links
        }

        self.business_logic.update_collection_version(version.version_id, body)

        self.crossref_provider.fetch_metadata.assert_not_called()
        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertIsNotNone(updated_version.publisher_metadata)

    def test_update_collection_change_doi(self):
        """
        A collection updated with a new DOI should get new publisher metadata from Crossref
        """
        metadata = self.sample_collection_metadata
        links = [Link("test doi", "DOI", "http://test.doi")]
        metadata.links = links

        self.crossref_provider.fetch_metadata = Mock(return_value={"authors": ["Test Author"]})

        version = self.initialize_unpublished_collection(metadata=metadata)
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

        self.business_logic.update_collection_version(version.version_id, body)

        self.crossref_provider.fetch_metadata.assert_called_once()
        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.publisher_metadata, expected_updated_publisher_metadata)


class TestUpdateCollectionDatasets(BaseBusinessLogicTestCase):

    def test_add_dataset_to_unpublished_collection_ok(self):
        """
        A dataset can be added to a collection when `ingest_dataset` is called.
        The resulting dataset should be empty and in a state ready for processing.
        """
        version = self.initialize_empty_unpublished_collection()
        url = "http://test/dataset.url"

        new_dataset_version_id = self.business_logic.ingest_dataset(version.version_id, url, None)

        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.processing_status, DatasetStatus.WAITING)

    def test_add_dataset_to_published_collection_fail(self):
        """
        Adding a dataset to a published collection should result in a failure
        """
        version = self.initialize_published_collection()
        url = "http://test/dataset.url"

        with self.assertRaises(DatasetIngestException):
            self.business_logic.ingest_dataset(version.version_id, url, None)

    def test_remove_dataset_from_unpublished_collection_ok(self):
        """
        A dataset can be removed from a collection version using `delete_dataset`.
        This should NOT delete the dataset but rather just update the collection_version -> dataset_version mapping
        """
        version = self.initialize_unpublished_collection()
        self.assertEqual(2, len(version.datasets))
        dataset_version_to_delete_id = version.datasets[0].version_id

        self.business_logic.delete_dataset(dataset_version_to_delete_id)

        new_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(0, len(new_version.datasets))

        # The dataset should still be present in the persistence store
        deleted_dataset = self.database_provider.get_dataset_version(dataset_version_to_delete_id)
        self.assertIsNotNone(deleted_dataset)

    def test_remove_dataset_from_published_collection_fail(self):
        """
        Removing a dataset from a published collection should fail
        """
        version = self.initialize_published_collection()
        self.assertEqual(2, len(version.datasets))
        dataset_version_to_delete_id = version.datasets[0].version_id

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.delete_dataset(dataset_version_to_delete_id)

    def test_replace_dataset_in_unpublished_collection_ok(self):
        """
        A dataset can be replaced from a collection version by calling `ingest_dataset`
        and specifying an existing dataset_version_id
        """
        version = self.initialize_unpublished_collection()
        dataset_version_to_replace_id = version.datasets[0].version_id
        url = "http://test/dataset.url"

        new_dataset_version_id = self.business_logic.ingest_dataset(
            version.version_id, 
            url, 
            dataset_version_to_replace_id
        )

        # Verify that the replaced dataset is in the right status
        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.processing_status, DatasetStatus.WAITING)

        # Verify that the old dataset is still existent
        old_dataset_version = self.database_provider.get_dataset_version(dataset_version_to_replace_id)
        self.assertIsNotNone(old_dataset_version)

    def test_replace_dataset_in_published_collection_fail(self):
        """
        Replacing a dataset that belongs to a published collection should fail
        """
        version = self.initialize_published_collection()
        dataset_version_to_replace_id = version.datasets[0].version_id
        url = "http://test/dataset.url"

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.ingest_dataset(
                version.version_id, 
                url, 
                dataset_version_to_replace_id
            )


class TestGetDataset(BaseBusinessLogicTestCase):

    def test_get_all_datasets_ok(self):
        """
        All dataset that belong to a published collection can be retrieved with `get_all_datasets`
        """
        # This will add 4 datasets, but only 2 should be retrieved by `get_all_datasets`
        published_version = self.initialize_published_collection() 
        unpublished_version = self.initialize_unpublished_collection()

        datasets = list(self.business_logic.get_all_datasets())
        self.assertEqual(2, len(datasets))
        self.assertCountEqual([d.version_id for d in datasets], [d.version_id for d in published_version.datasets])
        

    def test_get_dataset_artifacts_ok(self):
        """
        Artifacts belonging to a dataset can be obtained with `get_dataset_artifacts`
        """
        published_version = self.initialize_published_collection()
        dataset_id = published_version.datasets[0].dataset_id

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_id))
        self.assertEqual(3, len(artifacts))
        self.assertCountEqual([a.type for a in artifacts], ["H5AD", "CXG", "RDS"])

    def test_get_dataset_artifact_download_data_ok(self):
        """
        Calling `get_dataset_artifact_download_data` should yield downloadable data
        """
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        artifact = next(artifact for artifact in dataset.artifacts if artifact.type == "H5AD")
        self.assertIsNotNone(artifact)

        # TODO: requires mocking of the S3 provider. implement later
        download_data = self.business_logic.get_dataset_artifact_download_data(dataset.dataset_id, artifact.id)
        expected_download_data = DatasetArtifactDownloadData(
            "local.h5ad",
            "H5AD",
            0,
            ""
            )
        self.assertEqual(download_data, expected_download_data)


    def test_get_dataset_status_for_uploaded_dataset_ok(self):
        """
        Calling `get_dataset_status` should yield the processing status
        """
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        status = self.business_logic.get_dataset_status(dataset.dataset_id)
        self.assertEqual(status, DatasetStatus.UPLOADED)


class TestUpdateDataset(BaseBusinessLogicTestCase):

    def test_update_dataset_status_ok(self):
        """
        The dataset processing status can be updated using `update_dataset_status`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.business_logic.update_dataset_version_status(dataset.version_id, DatasetStatus.UPLOADED)
            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(version_from_db.processing_status, DatasetStatus.UPLOADED)

    def test_add_dataset_artifact_ok(self):
        """
        A dataset artifact can be added using `add_dataset_artifact`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.assertEqual(dataset.artifacts, [])
            self.business_logic.add_dataset_artifact(dataset.version_id, "H5AD", "http://fake.uri/artifact.h5ad")

            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(1, len(version_from_db.artifacts))
            self.assertEqual(version_from_db.artifacts[0].type, "H5AD")
            self.assertEqual(version_from_db.artifacts[0].uri, "http://fake.uri/artifact.h5ad")

        

    def test_add_dataset_artifact_wrong_type_fail(self):
        """
        Adding a dataset artifact with an unsupported type should fail
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.assertEqual(dataset.artifacts, [])
            with self.assertRaises(DatasetIngestException):
                self.business_logic.add_dataset_artifact(dataset.version_id, "BAD_TYPE", "http://fake.uri/artifact.h5ad")


class TestCollectionOperations(BaseBusinessLogicTestCase):

    def test_create_collection_version_ok(self):
        """
        A collection version can be created using `create_collection_version`
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        # The new version has a different version_id
        self.assertNotEqual(published_collection.version_id, new_version.version_id)

        # The new version has the same collection_id, collection metadata, and datasets
        self.assertEqual(published_collection.collection_id, new_version.collection_id)
        self.assertEqual(published_collection.metadata, new_version.metadata)
        self.assertEqual(published_collection.datasets, new_version.datasets)

        # The new version should not be published (aka Private)
        self.assertIsNone(new_version.published_at)

        # get_collection still retrieves the original version
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertNotEqual(version.version_id, new_version.version_id)


    def test_delete_collection_version_ok(self):
        """
        A collection version can be deleted using `delete_collection_version`
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        self.business_logic.delete_collection_version(new_version.version_id)

        # The version should no longer exist
        deleted_version = self.database_provider.get_collection_version(new_version.version_id)
        self.assertIsNone(deleted_version)

        # The original version should not be affected
        original_version = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertIsNotNone(original_version)
        self.assertEqual(published_collection, original_version)

        # get_collection still retrieves the original version
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertNotEqual(version.version_id, new_version.version_id)


    def test_publish_version_fails_on_published_collection(self):
        """
        `publish_collection_version` should fail if called on a collection version that is already published.
        """
        published_collection = self.initialize_published_collection()
        self.assertIsNotNone(published_collection.published_at)

        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(published_collection.version_id)

    def test_publish_version_ok(self):
        """
        A collection version can be published using `publish_collection`
        """
        unpublished_collection = self.initialize_unpublished_collection()
        self.business_logic.publish_collection_version(unpublished_collection.version_id)

        published_version = self.database_provider.get_collection_version(unpublished_collection.version_id)
        self.assertIsNotNone(published_version.published_at) # TODO: ideally, do a date assertion here (requires mocking)
        self.assertEqual(published_version.collection_id, unpublished_collection.version_id)

        # get_collection retrieves the new version
        version = self.business_logic.get_published_collection_version(unpublished_collection.collection_id)
        self.assertEqual(version.version_id, published_version.version_id)
        self.assertNotEqual(version.version_id, unpublished_collection.version_id)

    def test_publish_collection_with_no_datasets_fail(self):
        """
        Publishing a collection version with no datasets should fail
        """
        unpublished_collection = self.initialize_empty_unpublished_collection()

        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(unpublished_collection.version_id)

    def test_publish_version_with_updated_metadata_ok(self):
        """
        Publishing a collection with updated metadata should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        body = {
            "name": "new collection name",
            "description": "new collection description",
            "contact_name": "new contact name",
            "contact_email": "new_email@czi.com",
        }

        self.business_logic.update_collection_version(new_version.version_id, body)
        self.business_logic.publish_collection_version(new_version.version_id)

        new_version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertIsNone(new_version_from_db.published_at)
        self.assertEqual(new_version.collection_id, new_version_from_db.collection_id)
        self.assertEqual(new_version.version_id, new_version_from_db.version_id)
        self.assertEqual(new_version.datasets, new_version_from_db.datasets)

    def test_publish_version_with_removed_dataset_ok(self):
        """
        Publishing a version with a dataset removed should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        dataset_version_to_remove = new_version.datasets[0].version_id
        dataset_version_to_keep = new_version.datasets[1].version_id

        self.business_logic.delete_dataset(dataset_version_to_remove)

        # The new version should have only one dataset (before publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(1, len(version_from_db.datasets))

        # The original version should still have two datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))

        # Get collection retrieves the original version (with two datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertEqual(2, len(version.datasets))

        self.business_logic.publish_collection_version(new_version.version_id)

        # The new version should have only one dataset (after publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(1, len(version_from_db.datasets))

        # The old version should still have two datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))

        # Get collection retrieves the new version (with one datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(1, len(version.datasets))
        self.assertEqual(version.datasets[0].version_id, dataset_version_to_keep)

    def test_publish_version_with_added_dataset_ok(self):
        """
        Publishing a version with a dataset added should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        added_dataset_version_id = self.business_logic.ingest_dataset(new_version.version_id, "http://fake.url", None)
        self.complete_dataset_processing_with_success(added_dataset_version_id)

        # The new version should have three datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(3, len(version_from_db.datasets))

        # The original version should still have two datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))

        # Get collection retrieves the original version (with two datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertEqual(2, len(version.datasets))

        self.business_logic.publish_collection_version(new_version.version_id)

        # The new version should have three datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(3, len(version_from_db.datasets))

        # The old version should still have two datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))

        # Get collection retrieves the new version (with three datasets, including the new one)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(3, len(version.datasets))
        self.assertIn(added_dataset_version_id, [d.version_id for d in version.datasets])

    def test_publish_version_with_replaced_dataset_ok(self):
        """
        Publishing a version with a dataset replaced should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        # We will replace the first dataset
        dataset_id_to_replace = published_collection.datasets[0].version_id
        dataset_id_to_keep = published_collection.datasets[1].version_id

        replaced_dataset_version_id = self.business_logic.ingest_dataset(
            new_version.version_id, 
            "http://fake.url", 
            dataset_id_to_replace
        )

        self.complete_dataset_processing_with_success(replaced_dataset_version_id)

        # The new version should have the correct datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([replaced_dataset_version_id, dataset_id_to_keep], [d.version_id for d in version_from_db.datasets])

        # The original version should have the correct datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([dataset_id_to_replace, dataset_id_to_keep], [d.version_id for d in version_from_db.datasets])

        # Get collection retrieves the original version (with two datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertEqual(2, len(version.datasets))
        self.assertCountEqual([dataset_id_to_replace, dataset_id_to_keep], [d.version_id for d in version.datasets])

        self.business_logic.publish_collection_version(new_version.version_id)

        # The new version should have the correct datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([replaced_dataset_version_id, dataset_id_to_keep], [d.version_id for d in version_from_db.datasets])

        # The old version should have the correct datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([dataset_id_to_replace, dataset_id_to_keep], [d.version_id for d in version_from_db.datasets])

        # Get collection retrieves the new version (with two datasets, including the replaced one)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(3, len(version.datasets))
        self.assertCountEqual([replaced_dataset_version_id, dataset_id_to_keep], [d.version_id for d in version_from_db.datasets])


if __name__ == '__main__':
    unittest.main()

