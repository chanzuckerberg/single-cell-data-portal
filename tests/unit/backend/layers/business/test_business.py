import os
import unittest
from datetime import datetime
from unittest.mock import Mock, patch
from uuid import uuid4

from backend.layers.business.business import (
    BusinessLogic,
    CollectionMetadataUpdate,
    CollectionQueryFilter,
    DatasetArtifactDownloadData,
)
from backend.layers.business.exceptions import (
    CollectionCreationException,
    CollectionPublishException,
    CollectionUpdateException,
    CollectionVersionException,
    DatasetIngestException,
    DatasetNotFoundException,
)
from backend.layers.common.entities import (
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetArtifactType,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
    OntologyTermId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.persistence.persistence_mock import DatabaseProviderMock
from backend.layers.thirdparty.crossref_provider import (
    CrossrefDOINotFoundException,
    CrossrefException,
    CrossrefProviderInterface,
)
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProviderInterface
from tests.unit.backend.layers.fixtures import test_user_name

test_curator_name = "Test User"


class BaseBusinessLogicTestCase(unittest.TestCase):
    sample_collection_metadata: CollectionMetadata
    sample_dataset_metadata: DatasetMetadata

    test_user_name = "test_user_1"
    test_curator_name = "Test User"

    @classmethod
    def setUpClass(cls) -> None:
        cls.run_as_integration = os.environ.get("INTEGRATION_TEST", "false").lower() == "true"
        if cls.run_as_integration:
            database_uri = os.environ.get("DB_URI", "postgresql://postgres:secret@localhost")
            cls.database_provider = DatabaseProvider(database_uri=database_uri)
            cls.database_provider._drop_schema()

    def setUp(self) -> None:
        if self.run_as_integration:
            self.database_provider._create_schema()
        else:
            self.database_provider = DatabaseProviderMock()

        # Mock CorporaConfig
        # TODO: deduplicate with base_api
        def mock_config_fn(name):
            if name == "upload_max_file_size_gb":
                return 30

        mock_config = patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
        mock_config.start()

        # TODO: also deduplicate with base test
        from backend.layers.common import validation

        validation.valid_consortia = {
            "Consortia 1",
            "Consortia 2",
            "Consortia 3",
            "Consortia 4",
        }

        # By default these do nothing. They can be mocked by single test cases.
        self.crossref_provider = CrossrefProviderInterface()
        self.step_function_provider = StepFunctionProviderInterface()
        self.step_function_provider.start_step_function = Mock()
        self.s3_provider = S3ProviderInterface()
        self.uri_provider = UriProviderInterface()
        self.uri_provider.validate = Mock(return_value=True)  # By default, every link should be valid
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

        self.business_logic = BusinessLogic(
            database_provider=self.database_provider,
            crossref_provider=self.crossref_provider,
            step_function_provider=self.step_function_provider,
            s3_provider=self.s3_provider,
            uri_provider=self.uri_provider,
        )

        self.sample_collection_metadata = CollectionMetadata(
            "test collection 1",
            "description of test collection 1",
            "scientist",
            "scientist@czi.com",
            [],
            ["Consortia 1", "Consortia 2"],
        )

        self.sample_dataset_metadata = DatasetMetadata(
            name="test_dataset_name",
            organism=[OntologyTermId(label="test_organism_label", ontology_term_id="test_organism_term_id")],
            tissue=[OntologyTermId(label="test_tissue_label", ontology_term_id="test_tissue_term_id")],
            assay=[OntologyTermId(label="test_assay_label", ontology_term_id="test_assay_term_id")],
            disease=[OntologyTermId(label="test_disease_label", ontology_term_id="test_disease_term_id")],
            sex=[OntologyTermId(label="test_sex_label", ontology_term_id="test_sex_term_id")],
            self_reported_ethnicity=[
                OntologyTermId(
                    label="test_self_reported_ethnicity_label", ontology_term_id="test_self_reported_ethnicity_term_id"
                )
            ],
            development_stage=[
                OntologyTermId(label="test_development_stage_label", ontology_term_id="test_development_stage_term_id")
            ],
            cell_type=[OntologyTermId(label="test_cell_type_label", ontology_term_id="test_cell_type_term_id")],
            cell_count=10,
            schema_version="3.0.0",
            mean_genes_per_cell=0.5,
            batch_condition=["test_batch_1", "test_batch_2"],
            suspension_type=["test_suspension_type"],
            donor_id=["test_donor_1"],
            is_primary_data="BOTH",
            x_approximate_distribution="normal",
        )

    def tearDown(self):
        if self.run_as_integration:
            self.database_provider._drop_schema()

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.run_as_integration:
            cls.database_provider._engine.dispose()

    def initialize_empty_unpublished_collection(
        self, owner: str = test_user_name, curator_name: str = test_curator_name
    ) -> CollectionVersion:
        """
        Initializes an unpublished collection to be used for testing, with no datasets
        """
        version = self.database_provider.create_canonical_collection(
            owner,
            curator_name,
            self.sample_collection_metadata,
        )
        return version

    def initialize_unpublished_collection(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        complete_dataset_ingestion: bool = True,
    ) -> CollectionVersionWithDatasets:
        """
        Initializes an unpublished collection to be used for testing, with two datasets.
        By default also completes dataset ingestion (normally, a process that would be done asynchonously).
        Pass `complete_dataset_ingestion=False` if you want to initialize datasets only.
        """
        version = self.initialize_empty_unpublished_collection(owner, curator_name)
        for _ in range(2):
            dataset_version = self.database_provider.create_canonical_dataset(
                version.version_id,
            )
            self.database_provider.set_dataset_metadata(dataset_version.version_id, self.sample_dataset_metadata)
            self.database_provider.add_dataset_to_collection_version_mapping(
                version.version_id, dataset_version.version_id
            )
            if complete_dataset_ingestion:
                self.complete_dataset_processing_with_success(dataset_version.version_id)
        return self.database_provider.get_collection_version_with_datasets(version.version_id)

    def initialize_published_collection(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        published_at: datetime = None,
    ) -> CollectionVersionWithDatasets:
        """
        Initializes a published collection to be used for testing, with a single dataset
        """
        published_at = published_at or datetime.utcnow()
        version = self.initialize_unpublished_collection(owner, curator_name)
        self.database_provider.finalize_collection_version(version.collection_id, version.version_id, published_at)
        return self.database_provider.get_collection_version_with_datasets(version.version_id)

    def complete_dataset_processing_with_success(self, dataset_version_id: DatasetVersionId) -> None:
        """
        Test method that "completes" a dataset processing. This is necessary since dataset ingestion
        is a complex process which happens asynchronously, and cannot be easily mocked.
        """
        self.database_provider.add_dataset_artifact(
            dataset_version_id, DatasetArtifactType.H5AD.value, "s3://fake-bucket/local.h5ad"
        )
        self.database_provider.add_dataset_artifact(
            dataset_version_id, DatasetArtifactType.CXG.value, "s3://fake-bucket/local.cxg"
        )
        self.database_provider.add_dataset_artifact(
            dataset_version_id, DatasetArtifactType.RDS.value, "s3://fake-bucket/local.rds"
        )
        self.database_provider.update_dataset_upload_status(dataset_version_id, DatasetUploadStatus.UPLOADED)
        self.database_provider.update_dataset_validation_status(dataset_version_id, DatasetValidationStatus.VALID)
        self.database_provider.update_dataset_processing_status(dataset_version_id, DatasetProcessingStatus.SUCCESS)
        # TODO: if required, set the conversion status as well


class TestCreateCollection(BaseBusinessLogicTestCase):
    def test_create_collection_ok(self):
        """
        A collection can be created using `create_collection`
        """
        collection = self.business_logic.create_collection(
            test_user_name, test_curator_name, self.sample_collection_metadata
        )
        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(collection, collection_from_database)

    def test_create_collection_with_consortia_ok(self):
        """
        A collection with consortia can be created using `create_collection`
        """
        good_consortia = ["Consortia 1", "Consortia 2"]
        self.sample_collection_metadata.consortia = good_consortia
        collection = self.business_logic.create_collection(
            test_user_name, test_curator_name, self.sample_collection_metadata
        )
        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(good_consortia, collection_from_database.metadata.consortia)

    def test_create_collection_with_links_ok(self):
        """
        A collection with links can be created using `create_collection`
        """
        good_links = [
            Link("test link 1", "protocol", "http://example.com/protocol"),
            Link("test link 2", "other", "http://example.com/other"),
            Link(None, "other", "http://example.com/other"),  # names can be optional
        ]
        self.sample_collection_metadata.links = good_links
        collection = self.business_logic.create_collection(
            test_user_name, test_curator_name, self.sample_collection_metadata
        )
        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(good_links, collection_from_database.metadata.links)

    def test_create_collection_with_bad_links_fail(self):
        """
        A collection where a link isn't a http URI cannot be created
        """
        bad_links = [Link("test bad link", "other", "incorrect_url")]
        self.sample_collection_metadata.links = bad_links

        with self.assertRaises(CollectionCreationException) as ex:
            self.business_logic.create_collection(test_user_name, test_curator_name, self.sample_collection_metadata)

        self.assertEqual(
            ex.exception.errors, [{"name": "links[0]", "reason": "Invalid URL.", "value": "incorrect_url"}]
        )

    def test_create_collection_with_valid_doi_ok(self):
        """
        A collection can be created with a valid DOI, and the publisher metadata will be added
        to the collection
        """
        links_with_doi = [Link("test doi", "DOI", "http://good.doi")]
        self.sample_collection_metadata.links = links_with_doi

        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publisher_metadata)

        collection = self.business_logic.create_collection(
            test_user_name, test_curator_name, self.sample_collection_metadata
        )

        self.crossref_provider.fetch_metadata.assert_called_with("http://good.doi")

        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(1, len(collection_from_database.metadata.links))
        self.assertEqual(collection_from_database.metadata.links[0].uri, "http://good.doi")
        self.assertIsNotNone(collection_from_database.publisher_metadata)
        self.assertEqual(collection_from_database.publisher_metadata, expected_publisher_metadata)

    def test_create_collection_with_not_found_doi_fail(self):
        """
        A collection with a DOI that cannot be found on Crossref will not be created
        """
        links_with_doi = [Link("test doi", "DOI", "http://bad.doi")]

        self.sample_collection_metadata.links = links_with_doi

        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefDOINotFoundException("Error!"))

        with self.assertRaises(CollectionCreationException):
            self.business_logic.create_collection(test_user_name, test_curator_name, self.sample_collection_metadata)

    def test_create_collection_with_doi_error_ignores_metadata_ok(self):
        """
        A collection with an invalid DOI will be created with empty publisher metadata
        """
        links_with_doi = [Link("test doi", "DOI", "http://bad.doi")]

        self.sample_collection_metadata.links = links_with_doi

        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefException("Error!"))


class TestGetCollectionVersion(BaseBusinessLogicTestCase):
    def test_get_published_collection_version_for_published_collection_ok(self):
        """
        A published collection can be obtained using `get_collection`
        """
        version = self.initialize_published_collection()

        fetched_version = self.business_logic.get_published_collection_version(version.collection_id)

        self.assertIsNotNone(fetched_version)
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

        self.assertIsNone(fetched_version.published_at)
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
        All the collection versions should be returned by `get_collections`, including published, unpublished,
        and all owners
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
        If published filter flag is True, only published, non-tombstoned collections should be returned
        """
        self.initialize_unpublished_collection()
        self.initialize_published_collection()

        # Add a tombstoned Collection
        collection_version_to_tombstone: CollectionVersionWithDatasets = self.initialize_published_collection()
        self.database_provider.delete_canonical_collection(collection_version_to_tombstone.collection_id)

        # Confirm tombstoned Collection is in place
        all_collections_including_tombstones = self.database_provider.get_all_collections_versions(get_tombstoned=True)
        self.assertEqual(3, len(list(all_collections_including_tombstones)))

        all_collections_no_tombstones = self.database_provider.get_all_collections_versions()
        self.assertEqual(2, len(list(all_collections_no_tombstones)))

        filter = CollectionQueryFilter(is_published=True)
        versions = list(self.business_logic.get_collections(filter))

        self.assertEqual(1, len(versions))
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

        body = CollectionMetadataUpdate(
            name="new collection name",
            description="new collection description",
            contact_name="new contact name",
            contact_email="new_email@czi.com",
            links=[Link("test link 2", "other", "http://example.com/other")],
            consortia=["Consortia 1"],
        )

        self.business_logic.update_collection_version(version.version_id, body)

        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.metadata.name, body.name)
        self.assertEqual(updated_version.metadata.description, body.description)
        self.assertEqual(updated_version.metadata.contact_name, body.contact_name)
        self.assertEqual(updated_version.metadata.contact_email, body.contact_email)
        self.assertEqual(updated_version.metadata.links, body.links)
        self.assertEqual(updated_version.metadata.consortia, body.consortia)

    def test_update_collection_partial_ok(self):
        """
        `update_collection` should support partial updates: if only a subset of the fields are specified,
        the previous field should not be modified
        """
        version = self.initialize_unpublished_collection()

        body = CollectionMetadataUpdate(
            name="new collection name",
            description="new collection description",
            contact_name=None,
            contact_email=None,
            links=None,
            consortia=None,
        )

        self.business_logic.update_collection_version(version.version_id, body)

        updated_version = self.database_provider.get_collection_version(version.version_id)

        # These two fields should be updated
        self.assertEqual(updated_version.metadata.name, body.name)
        self.assertEqual(updated_version.metadata.description, body.description)

        # These three fields should remain the same
        self.assertEqual(updated_version.metadata.contact_name, self.sample_collection_metadata.contact_name)
        self.assertEqual(updated_version.metadata.contact_email, self.sample_collection_metadata.contact_email)
        self.assertEqual(updated_version.metadata.links, self.sample_collection_metadata.links)
        self.assertEqual(updated_version.metadata.consortia, self.sample_collection_metadata.consortia)

    def test_update_published_collection_fail(self):
        """
        Updating a collection version that is published should fail
        """
        version = self.initialize_published_collection()
        body = CollectionMetadataUpdate(
            name="new collection name",
            description="new collection description",
            contact_name="new contact name",
            contact_email="new_email@czi.com",
            links=None,
            consortia=None,
        )

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.update_collection_version(version.version_id, body)

    def test_update_collection_same_doi(self):
        """
        A collection updated with the same DOI should not trigger a Crossref call
        """
        metadata = self.sample_collection_metadata
        links = [Link("test doi", "DOI", "http://test.doi")]
        metadata.links = links

        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publisher_metadata)

        # We need to call `business_logic.create_collection` so that the publisher metadata is populated
        version = self.business_logic.create_collection(test_user_name, test_curator_name, metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = CollectionMetadataUpdate(
            name=None,
            description=None,
            contact_name=None,
            contact_email=None,
            links=links,
            consortia=None,
        )

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

        # We need to call `business_logic.create_collection` so that the publisher metadata is populated
        version = self.business_logic.create_collection(test_user_name, test_curator_name, metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = CollectionMetadataUpdate(
            name=None,
            description=None,
            contact_name=None,
            contact_email=None,
            links=[Link("new test doi", "DOI", "http://new.test.doi")],
            consortia=None,
        )

        expected_updated_publisher_metadata = {"authors": ["New Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_updated_publisher_metadata)

        self.business_logic.update_collection_version(version.version_id, body)

        self.crossref_provider.fetch_metadata.assert_called_once()
        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.publisher_metadata, expected_updated_publisher_metadata)


class TestUpdateCollectionDatasets(BaseBusinessLogicTestCase):
    def test_add_empty_dataset_ok(self):
        """
        An empty dataset can be added to a collection when `create_empty_dataset` is called.
        The resulting dataset should be empty and in a state ready for processing.
        """
        version = self.initialize_empty_unpublished_collection()

        new_dataset_version_id, _ = self.business_logic.create_empty_dataset(version.version_id)

        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, version.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.NA)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)

    def test_add_dataset_to_unpublished_collection_ok(self):
        """
        A dataset can be added to a collection when `ingest_dataset` is called.
        The resulting dataset should be empty and in a state ready for processing.
        """
        version = self.initialize_empty_unpublished_collection()
        url = "http://test/dataset.url"

        new_dataset_version_id, _ = self.business_logic.ingest_dataset(version.version_id, url, None, None)

        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, version.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.WAITING)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)

        self.step_function_provider.start_step_function.assert_called_once()

    def test_add_dataset_to_non_existing_collection_fail(self):
        """
        Calling `ingest_dataset` on a collection that does not exist should fail
        """
        url = "http://test/dataset.url"
        fake_collection_version_id = CollectionVersionId(str(uuid4()))

        with self.assertRaises(CollectionUpdateException) as ex:
            self.business_logic.ingest_dataset(fake_collection_version_id, url, None, None)
        self.assertEqual(ex.exception.errors, [f"Collection version {fake_collection_version_id} does not exist"])

        self.step_function_provider.start_step_function.assert_not_called()

    def test_add_dataset_to_published_collection_fail(self):
        """
        Adding a dataset to a published collection should result in a failure
        """
        version = self.initialize_published_collection()
        url = "http://test/dataset.url"

        with self.assertRaises(CollectionUpdateException) as ex:
            self.business_logic.ingest_dataset(version.version_id, url, None, None)
        self.assertEqual(ex.exception.errors, [f"Collection version {version.version_id.id} is published"])

        self.step_function_provider.start_step_function.assert_not_called()

    def test_add_dataset_with_invalid_link_fail(self):
        """
        Adding a dataset with a wrong link should result in a failure
        """
        version = self.initialize_published_collection()
        url = "http://bad.url"

        self.uri_provider.validate = Mock(return_value=False)

        with self.assertRaises(DatasetIngestException) as ex:
            self.business_logic.ingest_dataset(version.version_id, url, None, None)
        self.assertEqual(str(ex.exception), "Trying to upload invalid URI: http://bad.url")

        self.step_function_provider.start_step_function.assert_not_called()

    def test_remove_dataset_from_unpublished_collection_ok(self):
        """
        A dataset can be removed from a collection version using `delete_dataset`.
        This should NOT delete the dataset but rather just update the collection_version -> dataset_version mapping
        """
        version = self.initialize_unpublished_collection()
        self.assertEqual(2, len(version.datasets))
        dataset_version_to_delete_id = version.datasets[0].version_id

        self.business_logic.remove_dataset_version(version.version_id, dataset_version_to_delete_id)

        new_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(1, len(new_version.datasets))

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
            self.business_logic.remove_dataset_version(version.version_id, dataset_version_to_delete_id)

    def test_replace_dataset_in_unpublished_collection_ok(self):
        """
        A dataset can be replaced from a collection version by calling `ingest_dataset`
        and specifying an existing dataset_version_id
        """
        version = self.initialize_unpublished_collection()
        dataset_version_to_replace_id = version.datasets[0].version_id
        dataset_version_to_keep_id = version.datasets[1].version_id
        url = "http://test/dataset.url"

        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            version.version_id, url, None, dataset_version_to_replace_id
        )

        self.assertNotEqual(new_dataset_version_id, dataset_version_to_replace_id)

        # Verify that the replaced dataset is in the right status
        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.WAITING)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)

        # Verify that the old dataset is still existent
        old_dataset_version = self.database_provider.get_dataset_version(dataset_version_to_replace_id)
        self.assertIsNotNone(old_dataset_version)

        # Verify that the collection version points to the right datasets
        version_from_db = self.business_logic.get_collection_version(version.version_id)
        self.assertCountEqual(
            [d.version_id for d in version_from_db.datasets],
            [dataset_version_to_keep_id, new_dataset_version.version_id],
        )

        self.step_function_provider.start_step_function.assert_called_once()

    def test_replace_dataset_on_empty_dataset_does_not_generate_new_ok(self):
        """
        Calling `ingest_dataset` and specifying an existing dataset that is empty should start the ingestion
        without creating a new dataset version
        """
        version = self.initialize_empty_unpublished_collection()
        dataset_version_to_replace_id, _ = self.business_logic.create_empty_dataset(version.version_id)
        dataset_version_to_replace = self.business_logic.get_dataset_version(dataset_version_to_replace_id)
        url = "http://test/dataset.url"

        # Verify that the newly created dataset meets the condition for emptiness
        self.assertIsNotNone(dataset_version_to_replace)
        self.assertEqual(dataset_version_to_replace.status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.assertEqual(dataset_version_to_replace.artifacts, [])
        self.assertIsNone(dataset_version_to_replace.canonical_dataset.published_at)

        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            version.version_id, url, None, dataset_version_to_replace_id
        )

        self.assertEqual(new_dataset_version_id, dataset_version_to_replace_id)

        self.step_function_provider.start_step_function.assert_called_once()

    def test_replace_dataset_in_published_collection_fail(self):
        """
        Replacing a dataset that belongs to a published collection should fail
        """
        version = self.initialize_published_collection()
        dataset_version_to_replace_id = version.datasets[0].version_id
        url = "http://test/dataset.url"

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.ingest_dataset(version.version_id, url, None, dataset_version_to_replace_id)

        self.step_function_provider.start_step_function.assert_not_called()

    def test_replace_dataset_on_non_existing_dataset_fail(self):
        """
        Calling `ingest_dataset` and specifying a non existant dataset_version_id should fail
        """
        version = self.initialize_unpublished_collection()
        url = "http://test/dataset.url"

        with self.assertRaises(DatasetNotFoundException) as ex:
            self.business_logic.ingest_dataset(version.version_id, url, None, DatasetVersionId("fake_id"))
        self.assertEqual(str(ex.exception), "Dataset fake_id does not belong to the desired collection")

        self.step_function_provider.start_step_function.assert_not_called()

    def test_replace_dataset_in_wrong_status_fail(self):
        """
        Calling `ingest_dataset` and specifying an existing dataset should fail if the dataset is not in the
        final status
        """
        version = self.initialize_unpublished_collection()
        url = "http://test/dataset.url"

        dataset_version = version.datasets[0]
        self.database_provider.update_dataset_processing_status(
            dataset_version.version_id, DatasetProcessingStatus.PENDING
        )

        with self.assertRaises(DatasetIngestException) as ex:
            self.business_logic.ingest_dataset(version.version_id, url, None, dataset_version.version_id)
        self.assertEqual(
            str(ex.exception), f"Unable to reprocess dataset {dataset_version.version_id}: processing status is PENDING"
        )

        self.step_function_provider.start_step_function.assert_not_called()


class TestGetDataset(BaseBusinessLogicTestCase):
    def test_get_all_datasets_ok(self):
        """
        All dataset that belong to a published collection can be retrieved with `get_all_published_datasets`
        """
        # This will add 4 datasets, but only 2 should be retrieved by `get_all_datasets`
        published_version = self.initialize_published_collection()
        self.initialize_unpublished_collection()

        datasets = list(self.business_logic.get_all_published_datasets())
        self.assertEqual(2, len(datasets))
        self.assertCountEqual([d.version_id for d in datasets], [d.version_id for d in published_version.datasets])

    def test_get_dataset_artifacts_ok(self):
        """
        Artifacts belonging to a dataset can be obtained with `get_dataset_artifacts`
        """
        published_version = self.initialize_published_collection()
        dataset_version_id = published_version.datasets[0].version_id

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(3, len(artifacts))
        self.assertCountEqual(
            [a.type for a in artifacts], [DatasetArtifactType.H5AD, DatasetArtifactType.CXG, DatasetArtifactType.RDS]
        )
        self.assertCountEqual([a.get_file_name() for a in artifacts], ["local.h5ad", "local.cxg", "local.rds"])

    def test_get_dataset_artifact_download_data_ok(self):
        """
        Calling `get_dataset_artifact_download_data` should yield downloadable data
        """
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        artifact = next(artifact for artifact in dataset.artifacts if artifact.type == DatasetArtifactType.H5AD)
        self.assertIsNotNone(artifact)

        expected_file_size = 12345
        expected_presigned_url = "http://fake.presigned/url"

        self.s3_provider.get_file_size = Mock(return_value=expected_file_size)
        self.s3_provider.generate_presigned_url = Mock(return_value=expected_presigned_url)

        # TODO: requires mocking of the S3 provider. implement later
        download_data = self.business_logic.get_dataset_artifact_download_data(dataset.version_id, artifact.id)
        expected_download_data = DatasetArtifactDownloadData(
            "local.h5ad", DatasetArtifactType.H5AD, expected_file_size, expected_presigned_url
        )
        self.assertEqual(download_data, expected_download_data)

    def test_get_dataset_status_for_uploaded_dataset_ok(self):
        """
        Calling `get_dataset_status` should yield the processing status
        """
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        status = self.business_logic.get_dataset_status(dataset.version_id)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.SUCCESS)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
        self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)


class TestUpdateDataset(BaseBusinessLogicTestCase):
    def test_update_dataset_status_ok(self):
        """
        The dataset processing status can be updated using `update_dataset_status`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.business_logic.update_dataset_version_status(
                dataset.version_id, "upload", DatasetUploadStatus.UPLOADED
            )
            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(version_from_db.status.upload_status, DatasetUploadStatus.UPLOADED)

    def test_update_dataset_status_validation_message_ok(self):
        """
        The dataset validation message can be updated using `update_dataset_status`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.business_logic.update_dataset_version_status(
                dataset.version_id, "validation", DatasetValidationStatus.INVALID, "Validation error!"
            )
            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(version_from_db.status.validation_status, DatasetValidationStatus.INVALID)
            self.assertEqual(version_from_db.status.validation_message, "Validation error!")

    def test_add_dataset_artifact_ok(self):
        """
        A dataset artifact can be added using `add_dataset_artifact`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.assertEqual(dataset.artifacts, [])
            self.business_logic.add_dataset_artifact(dataset.version_id, "h5ad", "http://fake.uri/artifact.h5ad")

            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(1, len(version_from_db.artifacts))
            self.assertEqual(version_from_db.artifacts[0].type, DatasetArtifactType.H5AD.value)
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
                self.business_logic.add_dataset_artifact(
                    dataset.version_id, "BAD_TYPE", "http://fake.uri/artifact.h5ad"
                )

    def test_set_dataset_metadata_ok(self):
        """
        The dataset metadata can be set using `set_dataset_metadata`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.business_logic.set_dataset_metadata(dataset.version_id, self.sample_dataset_metadata)
            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(version_from_db.metadata, self.sample_dataset_metadata)


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

        # The canonical collection should be published
        self.assertIsNotNone(new_version.canonical_collection.originally_published_at)

        # get_collection still retrieves the original version
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertNotEqual(version.version_id, new_version.version_id)

    def test_create_collection_version_fails_if_other_versions(self):
        """
        A collection version can only be created if there are no other unpublished versions
        for that specific canonical collection
        """
        published_collection = self.initialize_published_collection()
        self.business_logic.create_collection_version(published_collection.collection_id)

        with self.assertRaises(CollectionVersionException):
            self.business_logic.create_collection_version(published_collection.collection_id)

    def test_create_collection_version_fails_if_collection_not_exists(self):
        """
        A collection version can only be created on an existing collection
        """
        non_existing_collection_id = CollectionId()
        with self.assertRaises(CollectionVersionException):
            self.business_logic.create_collection_version(CollectionId(non_existing_collection_id))

    def test_delete_collection_version_ok(self):
        """
        A collection version can be deleted using `delete_collection_version`
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        self.business_logic.delete_collection_version(new_version.version_id)

        # The version should no longer exist
        deleted_version = self.business_logic.get_collection_version(new_version.version_id)
        self.assertIsNone(deleted_version)

        # The original version should not be affected
        original_version = self.business_logic.get_collection_version(published_collection.version_id)
        self.assertIsNotNone(original_version)
        self.assertEqual(published_collection, original_version)

        # get_collection still retrieves the original version
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertNotEqual(version.version_id, new_version.version_id)

    def test_tombstone_collection_ok(self):
        """
        A collection can be marked as tombstoned using 'tombstone_collection'
        """
        published_collection = self.initialize_published_collection()
        self.business_logic.tombstone_collection(published_collection.collection_id)

        # The collection version canonical collection has tombstoned marked as True
        collection_version = self.business_logic.get_collection_version(published_collection.version_id)
        self.assertTrue(collection_version.canonical_collection.tombstoned)

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
        self.assertIsNotNone(
            published_version.published_at
        )  # TODO: ideally, do a date assertion here (requires mocking)
        self.assertIsNotNone(published_version.canonical_collection.originally_published_at)
        self.assertEqual(published_version.published_at, published_version.canonical_collection.originally_published_at)

        # The published and unpublished collection have the same collection_id and version_id
        self.assertEqual(published_version.collection_id, unpublished_collection.collection_id)
        self.assertEqual(published_version.version_id, unpublished_collection.version_id)

        # get_collection retrieves the correct version
        version = self.business_logic.get_published_collection_version(unpublished_collection.collection_id)
        if version:  # pylance
            self.assertEqual(version.version_id, published_version.version_id)

    def test_publish_collection_with_no_datasets_fail(self):
        """
        Publishing a collection version with no datasets should fail
        """
        unpublished_collection = self.initialize_empty_unpublished_collection()

        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(unpublished_collection.version_id)

    def test_publish_version_with_updated_metadata_ok(self):
        """
        Publishing a collection with updated metadata should succeed, but not update revised_at
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        body = CollectionMetadataUpdate(
            name="new collection name",
            description="new collection description",
            contact_name="new contact name",
            contact_email="new_email@czi.com",
            links=None,
            consortia=None,
        )

        self.assertIsNone(new_version.published_at)

        self.business_logic.update_collection_version(new_version.version_id, body)
        self.business_logic.publish_collection_version(new_version.version_id)

        new_version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        # Canonical Collection's revised_at should NOT be updated for collection metadata only revisions
        self.assertIsNone(new_version_from_db.canonical_collection.revised_at)
        self.assertEqual(published_collection.collection_id, new_version_from_db.collection_id)
        self.assertEqual(new_version.version_id, new_version_from_db.version_id)
        self.assertEqual([d.version_id for d in new_version.datasets], new_version_from_db.datasets)

        self.assertEqual(new_version_from_db.metadata.name, body.name)
        self.assertEqual(new_version_from_db.metadata.description, body.description)
        self.assertEqual(new_version_from_db.metadata.contact_name, body.contact_name)
        self.assertEqual(new_version_from_db.metadata.contact_email, body.contact_email)

    def test_publish_version_with_removed_dataset_ok(self):
        """
        Publishing a version with a dataset removed should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        dataset_version_to_remove = new_version.datasets[0]
        dataset_version_to_keep = new_version.datasets[1]

        self.business_logic.remove_dataset_version(new_version.version_id, dataset_version_to_remove.version_id)

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
        old_version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(old_version_from_db.datasets))

        # The Canonical Collection should have an updated revised_at timestamp
        self.assertEqual(version_from_db.canonical_collection.revised_at, version_from_db.published_at)
        self.assertTrue(version_from_db.canonical_collection.revised_at > old_version_from_db.published_at)

        # Get collection retrieves the new version (with one datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(1, len(version.datasets))
        self.assertEqual(version.datasets[0], dataset_version_to_keep)

    def test_publish_version_with_added_dataset_ok(self):
        """
        Publishing a version with a dataset added should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        added_dataset_version_id, _ = self.business_logic.ingest_dataset(
            new_version.version_id, "http://fake.url", None, None
        )
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
        old_version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(old_version_from_db.datasets))

        # The Canonical Collection should have an updated revised_at timestamp
        self.assertEqual(version_from_db.canonical_collection.revised_at, version_from_db.published_at)
        self.assertTrue(version_from_db.canonical_collection.revised_at > old_version_from_db.published_at)

        # Get collection retrieves the new version (with three datasets, including the new one)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(3, len(version.datasets))
        self.assertIn(added_dataset_version_id, [dataset.version_id for dataset in version.datasets])

    def test_publish_version_with_replaced_dataset_ok(self):
        """
        Publishing a version with a dataset replaced should succeed
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        # We will replace the first dataset
        dataset_id_to_replace = published_collection.datasets[0].version_id
        dataset_id_to_keep = published_collection.datasets[1].version_id

        replaced_dataset_version_id, _ = self.business_logic.ingest_dataset(
            new_version.version_id, "http://fake.url", None, dataset_id_to_replace
        )
        self.complete_dataset_processing_with_success(replaced_dataset_version_id)

        # The new version should have the correct datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([replaced_dataset_version_id, dataset_id_to_keep], version_from_db.datasets)

        # The original version should have the correct datasets (before publishing)
        version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([dataset_id_to_replace, dataset_id_to_keep], version_from_db.datasets)

        # Get collection retrieves the original version (with two datasets)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, published_collection.version_id)
        self.assertEqual(2, len(version.datasets))
        self.assertCountEqual(
            [dataset_id_to_replace, dataset_id_to_keep], [dataset.version_id for dataset in version.datasets]
        )

        self.business_logic.publish_collection_version(new_version.version_id)

        # The new version should have the correct datasets (after publishing)
        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        self.assertEqual(2, len(version_from_db.datasets))
        self.assertCountEqual([replaced_dataset_version_id, dataset_id_to_keep], version_from_db.datasets)

        # The old version should have the correct datasets (after publishing)
        old_version_from_db = self.database_provider.get_collection_version(published_collection.version_id)
        self.assertEqual(2, len(old_version_from_db.datasets))
        self.assertCountEqual([dataset_id_to_replace, dataset_id_to_keep], old_version_from_db.datasets)

        # The Canonical Collection should have an updated revised_at timestamp
        self.assertEqual(version_from_db.canonical_collection.revised_at, version_from_db.published_at)
        self.assertTrue(version_from_db.canonical_collection.revised_at > old_version_from_db.published_at)

        # Get collection retrieves the new version (with two datasets, including the replaced one)
        version = self.business_logic.get_published_collection_version(published_collection.collection_id)
        self.assertEqual(version.version_id, new_version.version_id)
        self.assertEqual(2, len(version.datasets))
        self.assertCountEqual(
            [replaced_dataset_version_id, dataset_id_to_keep], [dataset.version_id for dataset in version.datasets]
        )

    def test_publish_version_does_not_change_original_published_at_ok(self):
        """
        When publishing a collection version, the published_at for the canonical collection
        should not change
        """
        first_version = self.initialize_published_collection()
        second_version = self.business_logic.create_collection_version(first_version.collection_id)
        self.business_logic.publish_collection_version(second_version.version_id)

        canonical = self.business_logic.get_collection_version(second_version.version_id).canonical_collection
        self.assertEqual(canonical.originally_published_at, first_version.published_at)

    def test_get_all_collections_published_does_not_retrieve_old_versions(self):
        """
        `get_collections` with is_published=True should not return versions that were previously
        published but have a new version
        """
        first_version = self.initialize_published_collection()
        second_version = self.business_logic.create_collection_version(first_version.collection_id)
        self.business_logic.publish_collection_version(second_version.version_id)

        filter = CollectionQueryFilter(is_published=True)
        all_collections = list(self.business_logic.get_collections(filter))

        self.assertEqual(1, len(all_collections))
        self.assertEqual(all_collections[0].version_id, second_version.version_id)

        # The canonical collection published_at should point to the original publication time
        collection_version = self.business_logic.get_collection_version(all_collections[0].version_id)
        self.assertNotEqual(
            collection_version.canonical_collection.originally_published_at, all_collections[0].published_at
        )
        self.assertEqual(collection_version.canonical_collection.originally_published_at, first_version.published_at)

    def test_get_collection_versions_for_canonical_ok(self):
        """
        `get_collection_versions_from_canonical` can be used to retrieve all the versions
        for a canonical collection, whether those are published, active or neither
        """

        first_version = self.initialize_published_collection()
        second_version = self.business_logic.create_collection_version(first_version.collection_id)
        self.business_logic.publish_collection_version(second_version.version_id)
        third_version = self.business_logic.create_collection_version(first_version.collection_id)

        versions = self.business_logic.get_collection_versions_from_canonical(first_version.collection_id)
        versions = list(versions)

        self.assertEqual(3, len(versions))
        self.assertIn(first_version.version_id, [v.version_id for v in versions])
        self.assertIn(second_version.version_id, [v.version_id for v in versions])
        self.assertIn(third_version.version_id, [v.version_id for v in versions])

    def test_get_collection_version_from_canonical_published_ok(self):
        """
        `get_collection_version_from_canonical` retrieves the active published version connected
        to the canonical collection, when available
        """

        first_version = self.initialize_published_collection()
        self.business_logic.create_collection_version(first_version.collection_id)

        version = self.business_logic.get_collection_version_from_canonical(first_version.collection_id)
        self.assertIsNotNone(version)
        if version is not None:  # pylance
            self.assertEqual(version.version_id, first_version.version_id)
            self.assertIsNotNone(version.published_at)

    def test_get_collection_version_from_canonical_unpublished_ok(self):
        """
        `get_collection_version_from_canonical` retrieves the unpublished version connected
        to the canonical collection, if no published version is available
        """

        first_version = self.initialize_unpublished_collection()

        version = self.business_logic.get_collection_version_from_canonical(first_version.collection_id)
        self.assertIsNotNone(version)
        if version is not None:  # pylance
            self.assertEqual(version.version_id, first_version.version_id)
            self.assertIsNone(version.published_at)

    def test_dataset_published_at_with_new_dataset_ok(self):
        """
        When publishing a collection that contains a new dataset, `published_at` for the new dataset should be updated
        """

        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        added_dataset_version_id, _ = self.business_logic.ingest_dataset(
            new_version.version_id, "http://fake.url", None, None
        )
        self.complete_dataset_processing_with_success(added_dataset_version_id)

        self.business_logic.publish_collection_version(new_version.version_id)

        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        dataset_0 = self.database_provider.get_dataset_version(version_from_db.datasets[0])
        dataset_1 = self.database_provider.get_dataset_version(version_from_db.datasets[1])
        dataset_2 = self.database_provider.get_dataset_version(version_from_db.datasets[2])

        # dataset_2 is the new dataset
        # published_collection is the original, 2 dataset collection
        self.assertNotIn(dataset_2.version_id, [d.version_id for d in published_collection.datasets])
        self.assertIn(dataset_1.version_id, [d.version_id for d in published_collection.datasets])
        self.assertIn(dataset_0.version_id, [d.version_id for d in published_collection.datasets])

        # Published_at should be defined for all 3 datasets
        self.assertIsNotNone(dataset_0.canonical_dataset.published_at)
        self.assertIsNotNone(dataset_1.canonical_dataset.published_at)
        self.assertIsNotNone(dataset_2.canonical_dataset.published_at)

        self.assertGreater(dataset_2.canonical_dataset.published_at, published_collection.published_at)
        self.assertEqual(dataset_1.canonical_dataset.published_at, published_collection.published_at)
        self.assertEqual(dataset_0.canonical_dataset.published_at, published_collection.published_at)

    def test_dataset_published_at_with_replaced_dataset_ok(self):
        """
        When publishing a collection that contains a replaced dataset,
        `published_at` for the replaced dataset should not be updated
        """

        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        # We will replace the first dataset
        dataset_id_to_replace = published_collection.datasets[0].version_id
        dataset_id_to_keep = published_collection.datasets[1].version_id

        replaced_dataset_version_id, _ = self.business_logic.ingest_dataset(
            new_version.version_id, "http://fake.url", None, dataset_id_to_replace
        )

        self.complete_dataset_processing_with_success(replaced_dataset_version_id)

        self.business_logic.publish_collection_version(new_version.version_id)

        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        dataset_0 = self.database_provider.get_dataset_version(version_from_db.datasets[0])
        dataset_1 = self.database_provider.get_dataset_version(version_from_db.datasets[1])

        self.assertNotEqual(dataset_0.version_id, dataset_id_to_replace)
        self.assertEqual(dataset_1.version_id, dataset_id_to_keep)

        self.assertEqual(dataset_0.canonical_dataset.published_at, published_collection.published_at)
        self.assertEqual(dataset_1.canonical_dataset.published_at, published_collection.published_at)


if __name__ == "__main__":
    unittest.main()
