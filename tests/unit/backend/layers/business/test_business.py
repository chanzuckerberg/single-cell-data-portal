import itertools
import os
import unittest
import uuid
from copy import deepcopy
from datetime import datetime
from typing import List, Tuple
from unittest.mock import ANY, Mock, call
from uuid import uuid4

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from backend.common.corpora_config import CorporaConfig
from backend.common.providers.crossref_provider import (
    CrossrefDOINotFoundException,
    CrossrefException,
    CrossrefProviderInterface,
)
from backend.layers.business.business import (
    BusinessLogic,
    CollectionMetadataUpdate,
    CollectionQueryFilter,
    DatasetArtifactDownloadData,
)
from backend.layers.business.exceptions import (
    CollectionCreationException,
    CollectionDeleteException,
    CollectionIsPublishedException,
    CollectionPublishException,
    CollectionUpdateException,
    CollectionVersionException,
    DatasetIngestException,
    DatasetInWrongStatusException,
    DatasetIsTombstonedException,
    DatasetNotFoundException,
    InvalidIngestionManifestException,
    InvalidMetadataException,
    InvalidURIException,
    NoPreviousCollectionVersionException,
    NoPreviousDatasetVersionException,
)
from backend.layers.common.entities import (
    ARTIFACT_TO_EXTENSION,
    CollectionId,
    CollectionMetadata,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetId,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
    OntologyTermId,
    SpatialMetadata,
    TissueOntologyTermId,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.persistence.persistence_mock import DatabaseProviderMock
from backend.layers.thirdparty.batch_job_provider import BatchJobProviderInterface
from backend.layers.thirdparty.s3_exceptions import S3DeleteException
from backend.layers.thirdparty.s3_provider_mock import MockS3Provider
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
        # Set datasets bucket env var
        os.environ["DATASETS_BUCKET"] = "datasets"

    def setUp(self) -> None:
        if self.run_as_integration:
            self.database_provider._create_schema()
        else:
            self.database_provider = DatabaseProviderMock()

        # Mock CorporaConfig
        # TODO: deduplicate with base_api
        self.mock_config = CorporaConfig()
        self.mock_config.set(
            {
                "upload_max_file_size_gb": 30,
                "citation_update_feature_flag": "True",
                "dataset_assets_base_url": "https://dataset_assets_domain",
                "collections_base_url": "https://collections_domain",
            }
        )

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
        self.batch_job_provider = BatchJobProviderInterface()
        self.step_function_provider = StepFunctionProviderInterface()
        self.step_function_provider.start_step_function = Mock()
        self.s3_provider = MockS3Provider()
        self.uri_provider = UriProviderInterface()
        self.uri_provider.validate = Mock(return_value=True)  # By default, every link should be valid
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

        self.business_logic = BusinessLogic(
            database_provider=self.database_provider,
            batch_job_provider=self.batch_job_provider,
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
            tissue=[
                TissueOntologyTermId(
                    label="test_tissue_label", ontology_term_id="test_tissue_term_id", tissue_type="tissue"
                )
            ],
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
            primary_cell_count=5,
            schema_version="3.0.0",
            mean_genes_per_cell=0.5,
            batch_condition=["test_batch_1", "test_batch_2"],
            suspension_type=["test_suspension_type"],
            donor_id=["test_donor_1"],
            is_primary_data="BOTH",
            x_approximate_distribution="normal",
            default_embedding="X_embedding_1",
            embeddings=["X_embedding_1", "X_embedding_2"],
            feature_biotype=["gene"],
            feature_count=400,
            feature_reference=["NCBITaxon:9606"],
            raw_data_location="raw.X",
            citation="Publication: https://doi.org/12.2345/science.abc1234 Dataset Version: "
            "https://datasets.cellxgene.cziscience.com/dataset_id.h5ad curated and distributed by "
            "CZ CELLxGENE Discover in Collection: "
            "https://cellxgene.cziscience.com/collections/collection_id",
            spatial=SpatialMetadata(is_single=True, has_fullres=True),
        )
        self.s3_provider.mock_s3_fs = set()

    def tearDown(self):
        self.mock_config.reset()
        if self.run_as_integration:
            self.database_provider._drop_schema()

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.run_as_integration:
            cls.database_provider._engine.dispose()
        os.unsetenv("DATASETS_BUCKET")

    def initialize_empty_unpublished_collection(
        self, owner: str = test_user_name, curator_name: str = test_curator_name, links: List[Link] = None
    ) -> CollectionVersion:
        """
        Initializes an unpublished collection to be used for testing, with no datasets
        """
        if links:
            metadata = deepcopy(self.sample_collection_metadata)
            metadata.links = links
        else:
            metadata = self.sample_collection_metadata

        version = self.database_provider.create_canonical_collection(
            owner,
            curator_name,
            metadata,
        )
        return version

    def add_dataset_to_collection(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId = None,
        complete_dataset_ingestion: bool = True,
    ) -> None:
        """
        Links a dataset to a collection version
        """
        dataset_version_id = (
            self.database_provider.create_canonical_dataset(
                collection_version_id,
            ).version_id
            if dataset_version_id is None
            else dataset_version_id
        )
        self.database_provider.set_dataset_metadata(dataset_version_id, self.sample_dataset_metadata)
        self.database_provider.add_dataset_to_collection_version_mapping(collection_version_id, dataset_version_id)
        if complete_dataset_ingestion:
            self.complete_dataset_processing_with_success(dataset_version_id)
        return dataset_version_id

    def initialize_unpublished_collection(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        complete_dataset_ingestion: bool = True,
        num_datasets: int = 2,
        links: List[Link] = None,
    ) -> CollectionVersionWithDatasets:
        """
        Initializes an unpublished collection to be used for testing, with two datasets.
        By default also completes dataset ingestion (normally, a process that would be done asynchonously).
        Pass `complete_dataset_ingestion=False` if you want to initialize datasets only.
        """
        version = self.initialize_empty_unpublished_collection(owner, curator_name, links)
        for _ in range(num_datasets):
            self.add_dataset_to_collection(
                collection_version_id=version.version_id, complete_dataset_ingestion=complete_dataset_ingestion
            )
        return self.database_provider.get_collection_version_with_datasets(version.version_id)

    def initialize_published_collection(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        published_at: datetime = None,
        num_datasets: int = 2,
        links: List[Link] = None,
    ) -> CollectionVersionWithDatasets:
        """
        Initializes a published collection to be used for testing, with a single dataset
        """
        version = self.initialize_unpublished_collection(owner, curator_name, num_datasets=num_datasets, links=links)
        # published_at must be later than created_at for constituent Dataset Versions
        if published_at:
            assert published_at >= datetime.utcnow()
        else:
            published_at = datetime.utcnow()
        self.database_provider.finalize_collection_version(
            version.collection_id,
            version.version_id,
            "3.0.0",
            DATA_SUBMISSION_POLICY_VERSION,
            published_at=published_at,
        )
        return self.database_provider.get_collection_version_with_datasets(version.version_id)

    def initialize_collection_with_a_published_revision(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        published_at: datetime = None,
        num_datasets: int = 2,
    ) -> Tuple[CollectionVersionWithDatasets, CollectionVersionWithDatasets]:
        # Published with a published revision.
        published_version = self.initialize_published_collection(owner, curator_name, published_at, num_datasets)
        revision = self.business_logic.create_collection_version(published_version.collection_id)
        self.business_logic.publish_collection_version(revision.version_id)
        return (
            self.database_provider.get_collection_version_with_datasets(published_version.version_id),
            self.database_provider.get_collection_version_with_datasets(revision.version_id),
        )

    def initialize_collection_with_an_unpublished_revision(
        self,
        owner: str = test_user_name,
        curator_name: str = test_curator_name,
        published_at: datetime = None,
        num_datasets: int = 2,
        links: List[Link] = None,
    ) -> Tuple[CollectionVersionWithDatasets, CollectionVersionWithDatasets]:
        # Published with an unpublished revision
        published_version = self.initialize_published_collection(owner, curator_name, published_at, num_datasets, links)
        revision = self.business_logic.create_collection_version(published_version.collection_id)
        return (
            self.database_provider.get_collection_version_with_datasets(published_version.version_id),
            self.database_provider.get_collection_version_with_datasets(revision.version_id),
        )

    def complete_dataset_processing_with_success(self, dataset_version_id: DatasetVersionId, skip_atac=False) -> None:
        """
        Test method that "completes" a dataset processing. This is necessary since dataset ingestion
        is a complex process which happens asynchronously, and cannot be easily mocked.
        """

        def _add_artifact(bucket, key, key_type):
            ext = ARTIFACT_TO_EXTENSION[key_type]
            key_name = f"{key}.{ext}/" if key_type == DatasetArtifactType.CXG else f"{key}.{ext}"
            self.database_provider.create_dataset_artifact(dataset_version_id, key_type, f"s3://{bucket}/{key_name}")
            self.s3_provider.upload_file(None, bucket, key_name, None)

        _add_artifact("artifacts", f"{dataset_version_id}/raw", DatasetArtifactType.RAW_H5AD)
        # At present, not keeping public dataset assets as rows in DatasetArtifact table
        _add_artifact("datasets", f"{dataset_version_id}", DatasetArtifactType.H5AD)
        _add_artifact("datasets", f"{dataset_version_id}", DatasetArtifactType.RDS)
        _add_artifact("cellxgene", f"{dataset_version_id}", DatasetArtifactType.CXG)

        # special case for atac artifacts
        if not skip_atac:
            artifact_id = DatasetArtifactId()
            bucket = "datasets"

            key_name = f"{artifact_id}-fragment"
            _add_artifact(bucket, key_name, DatasetArtifactType.ATAC_FRAGMENT)
            _add_artifact(bucket, key_name, DatasetArtifactType.ATAC_INDEX)

        self.database_provider.update_dataset_upload_status(dataset_version_id, DatasetUploadStatus.UPLOADED)
        self.database_provider.update_dataset_validation_status(dataset_version_id, DatasetValidationStatus.VALID)
        self.database_provider.update_dataset_processing_status(dataset_version_id, DatasetProcessingStatus.SUCCESS)


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
        links_with_doi = [Link("test doi", "DOI", "https://doi.org/good/doi")]
        self.sample_collection_metadata.links = links_with_doi

        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, "good/doi", 17169328.664)
        )

        collection = self.business_logic.create_collection(
            test_user_name, test_curator_name, self.sample_collection_metadata
        )

        self.crossref_provider.fetch_metadata.assert_called_with("https://doi.org/good/doi")

        collection_from_database = self.database_provider.get_collection_version(collection.version_id)
        self.assertEqual(1, len(collection_from_database.metadata.links))
        self.assertEqual(collection_from_database.metadata.links[0].uri, "https://doi.org/good/doi")
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
        self.assertIsNotNone(fetched_version.data_submission_policy_version)

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
        self.assertIsNone(fetched_version.data_submission_policy_version)

    def test_get_collection_version_for_published_collection_ok(self):
        """
        A published collection version can be obtained using `get_collection_version`
        """
        version = self.initialize_published_collection()

        fetched_version = self.business_logic.get_collection_version(version.version_id)

        self.assertIsNotNone(fetched_version.published_at)
        self.assertEqual(fetched_version.metadata, version.metadata)

    def test_get_collection_version_for_tombstoned_collection(self):
        """
        A tombstoned Collection's versions can only be obtained if the `get_tombstoned` flag is explicitly set to True
        """
        version = self.initialize_published_collection()
        revision_version = self.business_logic.create_collection_version(version.collection_id)
        self.business_logic.publish_collection_version(revision_version.version_id)
        published_versions = self.business_logic.get_all_published_collection_versions_from_canonical(
            version.collection_id
        )
        self.assertEqual(2, len(list(published_versions)))

        # Tombstone the Collection
        self.business_logic.tombstone_collection(version.collection_id)

        with self.subTest("Does not retrieve tombstoned Collection version without get_tombstoned=True"):
            past_version_none = self.business_logic.get_collection_version(version.version_id)
            self.assertIsNone(past_version_none)

        with self.subTest("Retreives tombstoned Collection version when get_tombstoned=True"):
            past_version_tombstoned = self.business_logic.get_collection_version(
                version.version_id, get_tombstoned=True
            )
            self.assertIsNotNone(past_version_tombstoned)
            self.assertEqual(True, past_version_tombstoned.canonical_collection.tombstoned)

    def test_get_unpublished_collection_versions_from_canonical(self):
        """
        Unpublished Collection versions can be retrieved from a canonical Collection
        """
        published_collection = self.initialize_published_collection()
        non_migration_revision = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=False
        )
        migration_revision = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=True
        )

        unpublished_versions = self.business_logic.get_unpublished_collection_versions_from_canonical(
            published_collection.collection_id
        )
        unpublished_version_ids = [version.version_id.id for version in unpublished_versions]
        self.assertCountEqual(
            [non_migration_revision.version_id.id, migration_revision.version_id.id], unpublished_version_ids
        )


class TestGetAllCollections(BaseBusinessLogicTestCase):
    def test_get_all_collections_unfiltered_ok(self):
        """
        All the collection versions for all collections should be returned by `get_collections`, including published,
        unpublished, and all owners.
        """
        self.initialize_unpublished_collection()
        self.initialize_unpublished_collection(owner="test_user_2")
        self.initialize_published_collection(owner="test_user_2")
        self.initialize_collection_with_an_unpublished_revision()
        self.initialize_collection_with_a_published_revision()

        # TODO: this method should NOT be used without at least one filter. Maybe add an assertion to block it?
        filter = CollectionQueryFilter()
        versions = self.business_logic.get_collections(filter)

        self.assertEqual(7, len(list(versions)))

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
        self.initialize_collection_with_an_unpublished_revision()
        self.initialize_published_collection()
        self.initialize_collection_with_a_published_revision()

        # Add a tombstoned Collection
        collection_version_to_tombstone: CollectionVersionWithDatasets = self.initialize_published_collection()
        self.database_provider.tombstone_collection(collection_version_to_tombstone.collection_id)

        # Confirm tombstoned Collection is in place
        all_collections_including_tombstones = self.database_provider.get_all_collections_versions(get_tombstoned=True)
        self.assertEqual(7, len(list(all_collections_including_tombstones)))

        all_collections_no_tombstones = self.database_provider.get_all_collections_versions()
        self.assertEqual(6, len(list(all_collections_no_tombstones)))

        filter = CollectionQueryFilter(is_published=True)
        versions = list(self.business_logic.get_collections(filter))

        self.assertEqual(3, len(versions))
        for version in versions:
            self.assertIsNotNone(version.published_at)

    def test_get_all_collections_unpublished_ok(self):
        """
        If published filter flag is False, only unpublished collections should be returned
        """
        self.initialize_unpublished_collection()
        self.initialize_collection_with_an_unpublished_revision()
        self.initialize_published_collection()
        self.initialize_collection_with_a_published_revision()

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
        links = [Link("test doi", "DOI", "https://doi.org/test/doi")]
        metadata.links = links

        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, "test/doi", 17169328.664)
        )

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
        links = [Link("test doi", "DOI", "http://doi.org/test.doi")]
        metadata.links = links

        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, "test.doi", 17169328.664)
        )

        # We need to call `business_logic.create_collection` so that the publisher metadata is populated
        version = self.business_logic.create_collection(test_user_name, test_curator_name, metadata)
        self.crossref_provider.fetch_metadata.assert_called_once()
        self.crossref_provider.fetch_metadata.reset_mock()

        body = CollectionMetadataUpdate(
            name=None,
            description=None,
            contact_name=None,
            contact_email=None,
            links=[Link("new test doi", "DOI", "http://doi.org/new.test.doi")],
            consortia=None,
        )

        expected_updated_publisher_metadata = ({"authors": ["New Test Author"]}, "new.test.doi", 17169328.664)
        self.crossref_provider.fetch_metadata = Mock(return_value=expected_updated_publisher_metadata)
        self.batch_job_provider.start_metadata_update_batch_job = Mock()
        self.business_logic.update_collection_version(version.version_id, body)

        self.crossref_provider.fetch_metadata.assert_called_once()
        self.batch_job_provider.start_metadata_update_batch_job.assert_not_called()  # no datasets to update
        updated_version = self.database_provider.get_collection_version(version.version_id)
        self.assertEqual(updated_version.publisher_metadata, expected_updated_publisher_metadata[0])

    def test_update_collection_change_doi__trigger_dataset_artifact_updates(self):
        """
        A collection updated with a new DOI and containing datasets should trigger artifact updates to update citation
        for each dataset
        """
        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=2)

        body = CollectionMetadataUpdate(
            name=None,
            description=None,
            contact_name=None,
            contact_email=None,
            links=[Link("new test doi", "DOI", "http://doi.org/new.test.doi")],
            consortia=None,
        )
        self.batch_job_provider.start_metadata_update_batch_job = Mock()
        self.business_logic.generate_dataset_citation = Mock(return_value="test citation")
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["New Test Author"]}, "new.test.doi", 17169328.664)
        )

        self.business_logic.update_collection_version(revision.version_id, body)
        assert self.batch_job_provider.start_metadata_update_batch_job.call_count == 2
        self.batch_job_provider.start_metadata_update_batch_job.assert_has_calls(
            [
                call(revision.datasets[0].version_id, ANY, DatasetArtifactMetadataUpdate(citation="test citation")),
                call(revision.datasets[1].version_id, ANY, DatasetArtifactMetadataUpdate(citation="test citation")),
            ]
        )

    def test_update_published_collection_fail__dataset_status_pending(self):
        """
        Updating a collection version with a DOI change IF there is a dataset in non-finalized status should fail
        """
        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=2)
        self.business_logic.create_empty_dataset_version_for_current_dataset(
            revision.version_id,
            revision.datasets[0].version_id,
        )

        body = CollectionMetadataUpdate(
            name=None,
            description=None,
            contact_name=None,
            contact_email=None,
            links=[Link("new test doi", "DOI", "http://doi.org/new.test.doi")],
            consortia=None,
        )
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["New Test Author"]}, "new.test.doi", 17169328.664)
        )

        with self.assertRaises(CollectionUpdateException):
            self.business_logic.update_collection_version(revision.version_id, body)


class TestUpdateCollectionDatasets(BaseBusinessLogicTestCase):
    def test_add_empty_dataset_ok(self):
        """
        An empty dataset can be added to a collection when `create_empty_dataset` is called.
        The resulting dataset should be empty and in a state ready for processing.
        """
        version = self.initialize_empty_unpublished_collection()

        new_dataset_version_id = self.business_logic.create_empty_dataset(version.version_id).version_id

        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, version.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.NA)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)

    def test_create_empty_dataset_version_for_current_dataset(self):
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        new_dataset_version_id = self.business_logic.create_empty_dataset_version_for_current_dataset(
            revision.version_id, revision.datasets[0].version_id
        ).version_id
        revision_after_update = self.database_provider.get_collection_version(revision.version_id)
        self.assertEqual(len(revision_after_update.datasets), 1)
        self.assertEqual(revision_after_update.datasets[0], new_dataset_version_id)

        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, collection.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.WAITING)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)

    def test_create_empty_dataset_version_for_current_dataset__error_if_published(self):
        collection, _ = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        with self.assertRaises(CollectionIsPublishedException):
            self.business_logic.create_empty_dataset_version_for_current_dataset(
                collection.version_id, collection.datasets[0].version_id
            )

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

    def test_reingest_published_anndata_dataset(self):
        """A cellxgene public dataset url can be used to ingest a new dataset version."""

        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version = revision.datasets[0]
        url = f"https://dataset_assets_domain/{dataset_version.version_id}.h5ad"

        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, url, None, dataset_version.version_id
        )
        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, revision.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.WAITING)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.step_function_provider.start_step_function.assert_called_once_with(
            revision.version_id,
            new_dataset_version_id,
            f'{{"anndata":"s3://artifacts/{dataset_version.version_id}/raw.h5ad","atac_fragment":null}}',
        )

    def test_reingest_published_anndata_dataset__not_h5ad(self):
        """A cellxgene public dataset url used for reingesting an h5ad must be an h5ad"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version = revision.datasets[0]
        url = f"https://dataset_assets_domain/{dataset_version.version_id}.rds"
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, url, None, dataset_version.version_id)

    def test_reingest_published_anndata_dataset__not_in_found(self):
        """A cellxgene public dataset url must already be uploaded"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version_id = DatasetVersionId()
        url = f"https://dataset_assets_domain/{dataset_version_id.id}.h5ad"
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, url, None, dataset_version_id)

    def test_reginest_published_anndata_dataset__not_part_of_canonical_dataset(self):
        """A cellxgene public dataset url must be part of a version of the canonical dataset"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=2)
        dataset_version = revision.datasets[0]
        other_dataset_version = revision.datasets[1]
        url = f"https://dataset_assets_domain/{other_dataset_version.version_id}.h5ad"
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, url, None, dataset_version.version_id)

    def test_ingest_published_anndata_dataset_in_new_dataset__not_allowed(self):
        """A cellxgene public dataset url cannot be used to create a new canonical dataset."""
        published_dataset = self.initialize_published_collection().datasets[0]
        unpublished_collection = self.initialize_empty_unpublished_collection()
        url = f"https://dataset_assets_domain/{published_dataset.version_id.id}.h5ad"
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(unpublished_collection.version_id, url, None, None)

    def test_reingest_published_atac_dataset(self):
        """A cellxgene public dataset url can be used to ingest a new dataset version."""

        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version = revision.datasets[0]
        artifact_id = [a.id for a in dataset_version.artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
        anndata_url = f"https://dataset_assets_domain/{dataset_version.version_id}.h5ad"
        fragment_url = f"https://dataset_assets_domain/{artifact_id}-fragment.tsv.bgz"
        manifest = {"anndata": anndata_url, "atac_fragment": fragment_url}

        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, manifest, None, dataset_version.version_id
        )
        new_dataset_version = self.database_provider.get_dataset_version(new_dataset_version_id)
        self.assertIsNotNone(new_dataset_version)
        self.assertIsNone(new_dataset_version.metadata)
        self.assertEqual(new_dataset_version.collection_id, revision.collection_id)
        self.assertEqual(new_dataset_version.status.upload_status, DatasetUploadStatus.WAITING)
        self.assertEqual(new_dataset_version.status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.step_function_provider.start_step_function.assert_called_once_with(
            revision.version_id,
            new_dataset_version_id,
            f'{{"anndata":"s3://artifacts/{dataset_version.version_id}/raw.h5ad","atac_fragment":"{fragment_url}"}}',
        )

    def test_reingest_published_atac_dataset__not_atac(self):
        """A cellxgene public dataset url used for reingesting an h5ad must be an h5ad"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version = revision.datasets[0]
        anndata_url = f"https://dataset_assets_domain/{dataset_version.version_id}.h5ad"
        fragment_url = f"https://dataset_assets_domain/{dataset_version.version_id}.tsv"
        manifest = {"anndata": anndata_url, "atac_fragment": fragment_url}
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, manifest, None, dataset_version.version_id)

    def test_reingest_published_atac_dataset__not_in_found(self):
        """A cellxgene public dataset url must already be uploaded"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset_version = revision.datasets[0]
        anndata_url = f"https://dataset_assets_domain/{dataset_version.version_id}.h5ad"
        missing_dataset_version_id = DatasetVersionId()
        fragment_url = f"https://dataset_assets_domain/{revision.datasets[0].version_id}-fragment.tsv.bgz"
        manifest = {"anndata": anndata_url, "atac_fragment": fragment_url}
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, manifest, None, missing_dataset_version_id)

    def test_reginest_published_atac_dataset__not_part_of_canonical_dataset(self):
        """A cellxgene public dataset url must be part of a version of the canonical dataset"""
        collection, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=2)
        dataset_version = revision.datasets[0]
        anndata_url = f"https://dataset_assets_domain/{dataset_version.version_id}.h5ad"
        other_dataset_version = revision.datasets[1]
        fragment_url = f"https://dataset_assets_domain/{other_dataset_version.version_id}-fragment.tsv.bgz"
        manifest = {"anndata": anndata_url, "atac_fragment": fragment_url}
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(revision.version_id, manifest, None, dataset_version.version_id)

    def test_ingest_published_atac_dataset_in_new_dataset__not_allowed(self):
        """A cellxgene public dataset url cannot be used to create a new canonical dataset."""
        published_dataset = self.initialize_published_collection().datasets[0]
        unpublished_dataset = self.initialize_unpublished_collection().datasets[0]
        anndata_url = f"https://dataset_assets_domain/{unpublished_dataset.version_id}.h5ad"
        fragment_url = f"https://dataset_assets_domain/{published_dataset.version_id}-fragment.tsv.bgz"
        manifest = {"anndata": anndata_url, "atac_fragment": fragment_url}
        with self.assertRaises(InvalidIngestionManifestException):
            self.business_logic.ingest_dataset(unpublished_dataset.version_id, manifest, None, None)

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
        self.assertEqual(str(ex.exception), "Trying to upload invalid URI: http://bad.url/")

        self.step_function_provider.start_step_function.assert_not_called()

    def test_remove_dataset_from_unpublished_collection_ok(self):
        """
        A dataset can be removed from a collection version using `delete_dataset`. This should delete the dataset.
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

        # Verify that the previous dataset version is still present in database
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
        dataset_version_to_replace_id = self.business_logic.create_empty_dataset(version.version_id).version_id
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


class TestSetCollectionVersionDatasetsOrder(BaseBusinessLogicTestCase):
    def test_set_collection_version_datasets_order_ok(self):
        """
        The order of the datasets in a collection version is set using `set_collection_version_datasets_order`.
        """
        version = self.initialize_unpublished_collection(num_datasets=3)

        # Update cell counts to confirm custom order is returned and not default order.
        metadata_11 = deepcopy(self.sample_dataset_metadata)
        metadata_11.cell_count = 11
        self.database_provider.set_dataset_metadata(version.datasets[2].version_id, metadata_11)
        metadata_12 = deepcopy(self.sample_dataset_metadata)
        metadata_12.cell_count = 12
        self.database_provider.set_dataset_metadata(version.datasets[1].version_id, metadata_12)

        # Reverse and save the order of the dataset version IDs.
        dv_ids = [d.version_id for d in version.datasets]
        dv_ids.reverse()
        self.business_logic.set_collection_version_datasets_order(version.version_id, dv_ids)

        #  Confirm the saved order of collection version datasets in the database is correct.
        read_version = self.business_logic.get_collection_version(version.version_id)
        self.assertListEqual([dv.version_id for dv in read_version.datasets], dv_ids)

        # Confirm the collection version datasets are marked as custom ordered.
        self.assertTrue(read_version.has_custom_dataset_order)

    def test_set_collection_version_datasets_order_length_fail(self):
        """
        Attempting to set the order of the datasets in a collection version with a list of different length should fail.
        """
        version = self.initialize_unpublished_collection()

        # Remove a dataset version ID from the list to force length mismatch on save.
        dv_ids = [d.version_id for d in version.datasets]
        dv_ids.pop()

        with self.assertRaises(ValueError) as ex:
            self.business_logic.set_collection_version_datasets_order(version.version_id, dv_ids)
        self.assertEqual(
            str(ex.exception),
            f"Dataset Version IDs length does not match Collection Version {version.version_id.id} Datasets length",
        )

    def test_set_collection_version_datasets_order_invalid_dataset_fail(self):
        """
        Attempting to set the order of the datasets in a collection version with a list that contains an invalid
        dataset should fail.
        """
        version = self.initialize_unpublished_collection()

        # Update a dataset version ID in the list to force an invalid dataset on save.
        dv_ids = [d.version_id for d in version.datasets]
        dv_ids[0] = DatasetVersionId("fake_id")

        with self.assertRaises(ValueError) as ex:
            self.business_logic.set_collection_version_datasets_order(version.version_id, dv_ids)
        self.assertEqual(str(ex.exception), "Dataset Version IDs do not match saved Collection Version Dataset IDs")

    def test_set_collection_version_datasets_order_no_custom_order_ok(self):
        """
        The order of the datasets in a collection version is the default cell count, descending.
        """
        version = self.initialize_unpublished_collection(num_datasets=3)

        # Update cell counts to facilitate testing of order.
        metadata_11 = deepcopy(self.sample_dataset_metadata)
        metadata_11.cell_count = 11
        self.database_provider.set_dataset_metadata(version.datasets[2].version_id, metadata_11)
        metadata_12 = deepcopy(self.sample_dataset_metadata)
        metadata_12.cell_count = 12
        self.database_provider.set_dataset_metadata(version.datasets[1].version_id, metadata_12)

        # Confirm datasets on read collection version are ordered by cell count, descending.
        read_version = self.business_logic.get_collection_version(version.version_id)
        sorted_datasets = sorted(read_version.datasets, key=lambda d: d.metadata.cell_count, reverse=True)
        self.assertListEqual(read_version.datasets, sorted_datasets)

    def test_set_collection_version_datasets_order_replace_dataset_ok(self):
        """
        Replacing a dataset in a collection version does not alter the custom order of the datasets.
        """
        version = self.initialize_unpublished_collection(num_datasets=3)

        # Set custom order of dataset version IDs.
        dv_ids = [d.version_id for d in version.datasets]
        dv_ids.reverse()
        self.business_logic.set_collection_version_datasets_order(version.version_id, dv_ids)

        # Replace the first dataset in the collection version.
        dataset_id_to_replace = dv_ids[0]

        replaced_dataset_version_id, _ = self.business_logic.ingest_dataset(
            version.version_id, "http://fake.url", None, dataset_id_to_replace
        )
        self.business_logic.set_dataset_metadata(replaced_dataset_version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(replaced_dataset_version_id)

        # Create a list  of dataset version IDs with the replaced dataset version ID in the first position.
        updated_dv_ids = [replaced_dataset_version_id] + dv_ids[1:]

        # Confirm collection version lists datasets in the custom order.
        read_version = self.business_logic.get_collection_version(version.version_id)
        self.assertListEqual([dv.version_id for dv in read_version.datasets], updated_dv_ids)

        # Confirm the collection version datasets are marked as custom ordered.
        self.assertTrue(read_version.has_custom_dataset_order)


class TestDeleteDataset(BaseBusinessLogicTestCase):
    def test_delete_dataset_in_private_collection__ok(self):
        collection = self.initialize_unpublished_collection()
        dataset_version_id_to_remove = collection.datasets[0].version_id
        dataset_version_id_strs = [dv.version_id.id for dv in collection.datasets]
        self.assertIn(dataset_version_id_to_remove.id, dataset_version_id_strs)
        self.business_logic.remove_dataset_version(collection.version_id, dataset_version_id_to_remove)
        updated_collection = self.business_logic.get_collection_version(collection.version_id)
        dataset_version_id_strs = [dv.version_id.id for dv in updated_collection.datasets]
        self.assertNotIn(dataset_version_id_to_remove.id, dataset_version_id_strs)

    def test_delete_new_dataset_in_revision__ok(self):
        collection, revision = self.initialize_collection_with_an_unpublished_revision()
        new_dataset_version_id = self.business_logic.create_empty_dataset(revision.version_id).version_id
        revision = self.business_logic.get_collection_version(revision.version_id)
        self.complete_dataset_processing_with_success(new_dataset_version_id)
        # Dataset is added
        self.assertEqual(len(collection.datasets) + 1, len(revision.datasets))
        self.business_logic.remove_dataset_version(revision.version_id, new_dataset_version_id)
        revision = self.business_logic.get_collection_version(revision.version_id)
        # Dataset is removed
        self.assertEqual(len(collection.datasets), len(revision.datasets))

    def test_delete_published_dataset__fail(self):
        collection, revision = self.initialize_collection_with_an_unpublished_revision()
        with self.assertRaises(CollectionUpdateException):
            self.business_logic.remove_dataset_version(revision.version_id, revision.datasets[0].version_id)

    def test_delete_published_dataset_in_revision__success(self):
        collection, revision = self.initialize_collection_with_an_unpublished_revision()
        dataset_version_id_to_delete = revision.datasets[0].version_id
        self.business_logic.remove_dataset_version(
            revision.version_id, dataset_version_id_to_delete, delete_published=True
        )
        revision = self.business_logic.get_collection_version(revision.version_id)
        self.assertNotIn(dataset_version_id_to_delete.id, [dv.version_id.id for dv in revision.datasets])

    def test_cleanup_occurs_after_deletion_for_private_collection(self):
        """
        Updated Dataset assets and rows in private Collection should not be deleted until after Collection is deleted
        """
        collection = self.initialize_unpublished_collection(complete_dataset_ingestion=True)
        dataset_version_id_to_remove = collection.datasets[0].version_id
        self.business_logic.remove_dataset_version(collection.version_id, dataset_version_id_to_remove)
        self.assertIsNotNone(self.business_logic.get_dataset_version(dataset_version_id_to_remove))
        updated_collection = self.business_logic.get_collection_version(collection.version_id)
        dataset_version_id_strs = [dv.version_id.id for dv in updated_collection.datasets]
        self.assertNotIn(dataset_version_id_to_remove.id, dataset_version_id_strs)
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for d in collection.datasets for a in d.artifacts]
        self.business_logic.delete_collection_version(collection)
        self.assertIsNone(self.business_logic.get_dataset_version(dataset_version_id_to_remove))
        [self.assertFalse(self.s3_provider.uri_exists(a.uri)) for d in collection.datasets for a in d.artifacts]

    def test_cleanup_occurs_after_publication_of_private_collection(self):
        """
        Outdated, unused Dataset assets and rows in private Collection not deleted until after Collection is published
        """
        collection = self.initialize_unpublished_collection(complete_dataset_ingestion=True, num_datasets=3)
        # This dataset will be kept
        dataset_to_keep = collection.datasets[1]

        # Thisdataset will be replaced
        dataset_to_replace = collection.datasets[0]
        new_dataset_version = self.database_provider.replace_dataset_in_collection_version(
            collection.version_id, dataset_to_replace.version_id
        )
        self.business_logic.set_dataset_metadata(new_dataset_version.version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(new_dataset_version.version_id)
        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version.version_id)

        # This dataset will be checked for artifact retention across dataset versions
        dataset_to_check_artifact_retention = collection.datasets[2]
        new_atac_dataset_version = self.database_provider.replace_dataset_in_collection_version(
            collection.version_id, dataset_to_check_artifact_retention.version_id
        )
        self.business_logic.set_dataset_metadata(new_atac_dataset_version.version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(new_atac_dataset_version.version_id, skip_atac=True)
        # add fragment and index artifacts to the new_ATAC_dataset from the dataset_to_check_artifact_retention
        fragment_id = [
            a for a in dataset_to_check_artifact_retention.artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT
        ][0].id
        self.database_provider.add_artifact_to_dataset_version(new_atac_dataset_version.version_id, fragment_id)
        fragment_index_id = [
            a for a in dataset_to_check_artifact_retention.artifacts if a.type == DatasetArtifactType.ATAC_INDEX
        ][0].id
        self.database_provider.add_artifact_to_dataset_version(new_atac_dataset_version.version_id, fragment_index_id)
        new_atac_dataset_version = self.business_logic.get_dataset_version(new_atac_dataset_version.version_id)

        # Get Updates collection
        updated_collection = self.business_logic.get_collection_version(collection.version_id)

        # Verify that the old datasets have been replaced in the collection
        dataset_version_id_strs = [dv.version_id.id for dv in updated_collection.datasets]
        self.assertNotIn(dataset_to_replace.version_id.id, dataset_version_id_strs)
        self.assertNotIn(dataset_to_check_artifact_retention.version_id.id, dataset_version_id_strs)

        # Verify that the the old datasets are still present in the database
        self.assertIsNotNone(self.business_logic.get_dataset_version(dataset_to_replace.version_id))
        self.assertIsNotNone(self.business_logic.get_dataset_version(dataset_to_check_artifact_retention.version_id))

        # Verify all artifacts for the datasets exist before publication
        for artifact in itertools.chain(
            dataset_to_replace.artifacts,
            dataset_to_keep.artifacts,
            dataset_to_check_artifact_retention.artifacts,
            new_dataset_version.artifacts,
            new_atac_dataset_version.artifacts,
        ):
            self.assertTrue(self.s3_provider.uri_exists(artifact.uri))

        # Publish the collection
        self.business_logic.publish_collection_version(collection.version_id)

        # After publication, these dataset should be gone
        self.assertIsNone(self.business_logic.get_dataset_version(dataset_to_replace.version_id))
        self.assertIsNone(self.business_logic.get_dataset_version(dataset_to_check_artifact_retention.version_id))
        # Artifacts for dataset_to_replace should be gone
        [self.assertFalse(self.s3_provider.uri_exists(a.uri), f"Found {a.uri}") for a in dataset_to_replace.artifacts]
        # Some Artifact for dataset_to_check_artifact_retention should be gone
        retained_artifacts = [DatasetArtifactType.ATAC_FRAGMENT, DatasetArtifactType.ATAC_INDEX]
        [
            self.assertFalse(self.s3_provider.uri_exists(a.uri), f"Found {a.uri}")
            for a in dataset_to_check_artifact_retention.artifacts
            if a.type not in retained_artifacts
        ]
        [
            self.assertTrue(self.s3_provider.uri_exists(a.uri), f"Missing {a.uri}")
            for a in dataset_to_check_artifact_retention.artifacts
            if a.type in retained_artifacts
        ]
        # Artifacts for dataset_to_keep should remain
        [self.assertTrue(self.s3_provider.uri_exists(a.uri), f"Missing {a.uri}") for a in dataset_to_keep.artifacts]
        # Artifacts for new_dataset_version should remain
        [self.assertTrue(self.s3_provider.uri_exists(a.uri), f"Missing {a.uri}") for a in new_dataset_version.artifacts]
        # Artifacts for new_atac_dataset_version should remain
        [
            self.assertTrue(self.s3_provider.uri_exists(a.uri), f"Missing {a.uri}")
            for a in new_atac_dataset_version.artifacts
        ]

    def test_cleanup_occurs_after_revision_is_canceled(self):
        """
        New Dataset assets and rows not deleted until after revision is canceled (for previously-published Collection)
        """
        collection = self.initialize_published_collection()
        dataset_to_replace = collection.datasets[0]

        revision = self.business_logic.create_collection_version(collection.collection_id)

        updated_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, "http://fake.url", None, dataset_to_replace.version_id
        )
        self.business_logic.set_dataset_metadata(updated_dataset_version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(updated_dataset_version_id)
        updated_dataset_version = self.business_logic.get_dataset_version(updated_dataset_version_id)

        updated_updated_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, "http://fake.url", None, updated_dataset_version_id
        )
        self.business_logic.set_dataset_metadata(updated_updated_dataset_version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(updated_updated_dataset_version_id)
        updated_updated_dataset_version = self.business_logic.get_dataset_version(updated_updated_dataset_version_id)

        # Add new Dataset
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, "http://fake.url", None, None
        )
        self.business_logic.set_dataset_metadata(new_dataset_version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(new_dataset_version_id)
        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)

        # Update the new Dataset
        updated_new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, "http://fake.url", None, new_dataset_version_id
        )
        self.business_logic.set_dataset_metadata(updated_new_dataset_version_id, self.sample_dataset_metadata)
        self.complete_dataset_processing_with_success(updated_new_dataset_version_id)
        updated_new_dataset_version = self.business_logic.get_dataset_version(updated_new_dataset_version_id)

        # Artifacts for initial published Datasets should exist
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for dv in collection.datasets for a in dv.artifacts]
        # Artifacts for updated_dataset_version and updated_updated_dataset_version should exist
        [
            self.assertTrue(self.s3_provider.uri_exists(a.uri))
            for d in (updated_dataset_version, updated_updated_dataset_version)
            for a in d.artifacts
        ]
        # Artifacts for new_dataset_version should exist
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for a in new_dataset_version.artifacts]
        # Artifacts for updated_new_dataset_version should exist
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for a in updated_new_dataset_version.artifacts]

        self.business_logic.publish_collection_version(revision.version_id)

        # Initial published Datasets should remain
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for dv in collection.datasets for a in dv.artifacts]
        [self.assertIsNotNone(self.business_logic.get_dataset_version(dv.version_id)) for dv in collection.datasets]

        # updated_dataset_version should be gone
        [self.assertFalse(self.s3_provider.uri_exists(a.uri)) for a in updated_dataset_version.artifacts]
        self.assertIsNone(self.business_logic.get_dataset_version(updated_dataset_version_id))

        # updated_updated_dataset_version should remain
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for a in updated_updated_dataset_version.artifacts]
        self.assertIsNotNone(self.business_logic.get_dataset_version(updated_updated_dataset_version_id))

        # Artifacts for new_dataset_version should not exist
        [self.assertFalse(self.s3_provider.uri_exists(a.uri)) for a in new_dataset_version.artifacts]
        # Artifacts for updated_new_dataset_version should exist
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for a in updated_new_dataset_version.artifacts]

        self.assertIsNone(self.business_logic.get_dataset_version(new_dataset_version_id))
        self.assertIsNotNone(self.business_logic.get_dataset_version(updated_new_dataset_version_id))

    def test_delete_artifacts_delete_no_artifacts(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts
        s3_provider = self.business_logic.s3_provider

        # Arrange
        # Act
        delete_artifacts([])
        # Assert
        s3_provider.delete_prefix.assert_not_called()
        s3_provider.delete_files.assert_not_called()

    def test_delete_artifacts_artifact_type_doesnt_match_regex(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts
        s3_provider = self.business_logic.s3_provider
        bucket = "bucket"

        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.CXG, uri=f"s3://{bucket}/file.h5ad"),
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.H5AD, uri=f"s3://{bucket}/file.cxg/"),
        ]
        # Act
        delete_artifacts(artifacts)
        # Assert
        s3_provider.delete_prefix.assert_not_called()
        s3_provider.delete_files.assert_not_called()

    def test_delete_artifacts_delete_cxg(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts
        s3_provider = self.business_logic.s3_provider
        bucket = "bucket"

        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.CXG, uri=f"s3://{bucket}/file.cxg/"),
        ]
        # Act
        delete_artifacts(artifacts)
        # Assert
        s3_provider.delete_prefix.assert_called_with(bucket, "file.cxg/")

    def test_delete_artifacts_with_rdev_prefix(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts
        s3_provider = self.business_logic.s3_provider
        bucket = "bucket"

        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.CXG, uri=f"s3://{bucket}/rdev/file.cxg/"),
        ]
        # Act
        delete_artifacts(artifacts)
        # Assert
        s3_provider.delete_prefix.assert_called_with(bucket, "rdev/file.cxg/")

    def test_delete_artifacts_with_no_match(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts

        ## not matching regex
        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.H5AD, uri="s3://file.h5ad"),
        ]
        # Act + Assert
        self.assertRaises(InvalidURIException, delete_artifacts, artifacts)

    def test_delete_artifacts_delete_h5ad(self):
        self.business_logic.s3_provider = Mock()
        delete_artifacts = self.business_logic.delete_artifacts
        s3_provider = self.business_logic.s3_provider
        bucket = "bucket"

        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.H5AD, uri=f"s3://{bucket}/file.h5ad"),
        ]
        # Act
        delete_artifacts(artifacts)
        # Assert
        s3_provider.delete_files.assert_called_with(bucket, ["file.h5ad"])

    def test_delete_artifacts_CollectionDeleteException_with_H5AD(self):
        self.business_logic.s3_provider.delete_files = Mock(side_effect=S3DeleteException("error"))
        delete_artifacts = self.business_logic.delete_artifacts
        bucket = "bucket"

        # Arrange
        artifacts: List[DatasetArtifact] = [
            DatasetArtifact(id=DatasetArtifactId(), type=DatasetArtifactType.H5AD, uri=f"s3://{bucket}/file.h5ad"),
        ]

        # Act + Assert
        self.assertRaises(CollectionDeleteException, delete_artifacts, artifacts)

    def test_delete_dataset_versions(self):
        collection = self.initialize_unpublished_collection()
        dataset_to_delete = collection.datasets[0]
        uris_to_delete = [a.uri for a in dataset_to_delete.artifacts]
        dataset_artifact_ids = [a.id for a in dataset_to_delete.artifacts]

        self.business_logic.delete_dataset_versions([dataset_to_delete])

        self.assertIsNone(self.business_logic.get_dataset_version(dataset_to_delete.version_id))
        self.assertEqual(self.database_provider.get_dataset_artifacts(dataset_artifact_ids), [])
        [self.assertFalse(self.s3_provider.uri_exists(uri)) for uri in uris_to_delete]


class TestGetDataset(BaseBusinessLogicTestCase):
    def test_get_all_datasets_ok(self):
        """
        All dataset that belong to a published collection can be retrieved with `get_all_mapped_datasets`
        """
        # This will add 4 datasets, but only 2 should be retrieved by `get_all_datasets`
        published_version = self.initialize_published_collection()
        self.initialize_unpublished_collection()

        datasets = self.business_logic.get_all_mapped_datasets()
        self.assertEqual(2, len(datasets))
        self.assertCountEqual([d.version_id for d in datasets], [d.version_id for d in published_version.datasets])

    def test_get_dataset_version_from_canonical(self):
        """
        Get currently published dataset version using canonical dataset ID, or most recent active unpublished dataset
        if none are published.
        """
        with self.subTest("Dataset is published with a revision open, get published dataset version"):
            published_version = self.initialize_published_collection()
            published_dataset = self.business_logic.get_all_mapped_datasets()[0]
            self.business_logic.create_collection_version(published_version.collection_id)

            dataset_version = self.business_logic.get_dataset_version_from_canonical(published_dataset.dataset_id)
            self.assertEqual(dataset_version.version_id, published_dataset.version_id)
        with self.subTest("Dataset has never been published, get latest unpublished version"):
            unpublished_collection = self.initialize_unpublished_collection()
            init_dataset = unpublished_collection.datasets[0]
            new_dataset = self.database_provider.replace_dataset_in_collection_version(
                unpublished_collection.version_id, init_dataset.version_id
            )

            dataset_version = self.business_logic.get_dataset_version_from_canonical(init_dataset.dataset_id)
            self.assertEqual(dataset_version.version_id, new_dataset.version_id)
        with self.subTest("Dataset does not exist"):
            dataset_version = self.business_logic.get_dataset_version_from_canonical(DatasetId(str(uuid.uuid4())))
            self.assertIsNone(dataset_version)
        with self.subTest("DatasetVersion has been rolled back from most recently-created"):
            unpublished_collection = self.initialize_unpublished_collection()
            init_dataset = unpublished_collection.datasets[0]
            new_dataset = self.database_provider.replace_dataset_in_collection_version(
                unpublished_collection.version_id, init_dataset.version_id
            )
            self.business_logic.restore_previous_dataset_version(
                unpublished_collection.version_id, new_dataset.dataset_id
            )
            dataset_version = self.business_logic.get_dataset_version_from_canonical(init_dataset.dataset_id)
            self.assertEqual(dataset_version.version_id, init_dataset.version_id)

    def test_get_prior_published_versions_for_dataset(self):
        """
        Given a canonical dataset id, return all its DatasetVersions that have been part of published CollectionVersions
        """
        self.initialize_published_collection()
        published_version = self.initialize_published_collection()
        collection_id = published_version.collection_id
        dataset = published_version.datasets[0]
        # Revision 1 (to publish)
        revision_id = self.business_logic.create_collection_version(collection_id).version_id
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision_id, "http://fake.url", None, dataset.version_id
        )
        self.business_logic.set_dataset_metadata(new_dataset_version_id, self.sample_dataset_metadata)
        self.business_logic.publish_collection_version(revision_id)
        revision = self.business_logic.get_collection_version(revision_id)
        # Revision 2 (not to publish)
        unpublished_collection_version_id = self.business_logic.create_collection_version(collection_id).version_id
        unpublished_dataset_version_id, _ = self.business_logic.ingest_dataset(
            unpublished_collection_version_id, "http://fake.url", None, new_dataset_version_id
        )

        version_history = self.business_logic.get_prior_published_versions_for_dataset(dataset.dataset_id)
        # Check that only published datasets appear
        self.assertEqual(
            [dataset.version_id, new_dataset_version_id], [version.version_id for version in version_history]
        )
        self.assertEqual(
            [published_version.published_at, revision.published_at],
            [version.published_at for version in version_history],
        )

    def test_get_prior_published_dataset_version(self):
        """
        Given a dataset version id, return the DatasetVersion IF its been part of a published CollectionVersion
        """
        unpublished_dataset_version_id = self.initialize_unpublished_collection().datasets[0].version_id
        initial_published_collection_version = self.initialize_published_collection()
        collection_id = initial_published_collection_version.collection_id
        dataset = initial_published_collection_version.datasets[0]
        # Revision 1 (to publish)
        collection_revision_id = self.business_logic.create_collection_version(collection_id).version_id
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_revision_id, "http://fake.url", None, dataset.version_id
        )
        self.business_logic.set_dataset_metadata(new_dataset_version_id, self.sample_dataset_metadata)
        self.business_logic.publish_collection_version(collection_revision_id)
        # Revision 2 (not to publish)
        unpublished_collection_version_id = self.business_logic.create_collection_version(collection_id).version_id
        unpublished_dataset_revision_id, _ = self.business_logic.ingest_dataset(
            unpublished_collection_version_id, "http://fake.url", None, new_dataset_version_id
        )

        # test get unpublished dataset version
        unpublished_dataset_version = self.business_logic.get_prior_published_dataset_version(
            unpublished_dataset_version_id
        )
        self.assertIsNone(unpublished_dataset_version)
        # test get past published dataset version
        initial_published_dataset_version = self.business_logic.get_prior_published_dataset_version(dataset.version_id)
        self.assertEqual(initial_published_dataset_version.version_id, dataset.version_id)
        self.assertEqual(
            initial_published_dataset_version.collection_version_id, initial_published_collection_version.version_id
        )
        # test currently published dataset version
        published_dataset_revision = self.business_logic.get_prior_published_dataset_version(new_dataset_version_id)
        self.assertEqual(published_dataset_revision.version_id, new_dataset_version_id)
        self.assertEqual(published_dataset_revision.collection_version_id, collection_revision_id)
        # test get unpublished dataset revision
        unpublished_dataset_revision = self.business_logic.get_prior_published_dataset_version(
            unpublished_dataset_revision_id
        )
        self.assertIsNone(unpublished_dataset_revision)

    def test_get_dataset_version(self):
        collection = self.initialize_published_collection()
        dataset_to_tombstone = collection.datasets[0]

        revision = self.business_logic.create_collection_version(collection.collection_id)
        # Tombstone Dataset
        self.business_logic.remove_dataset_version(revision.version_id, dataset_to_tombstone.version_id, True)
        self.business_logic.publish_collection_version(revision.version_id)
        self.assertIsNone(self.business_logic.get_dataset_version(dataset_to_tombstone.version_id))
        self.assertIsNotNone(
            self.database_provider.get_dataset_version(dataset_to_tombstone.version_id, get_tombstoned=True)
        )

    def test_get_dataset_artifacts_ok(self):
        """
        Artifacts belonging to a dataset can be obtained with `get_dataset_artifacts`
        """
        published_version = self.initialize_published_collection()
        dataset_version_id = published_version.datasets[0].version_id

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(len(published_version.datasets[0].artifacts), len(artifacts))

    def test_get_dataset_artifact_download_data_ok(self):
        """
        Calling `get_dataset_artifact_download_data` should yield downloadable data
        """
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        artifact = next(artifact for artifact in dataset.artifacts if artifact.type == DatasetArtifactType.H5AD)
        self.assertIsNotNone(artifact)

        expected_file_size = 12345
        expected_permanent_url = "http://fake.permanent/url"

        self.s3_provider.get_file_size = Mock(return_value=expected_file_size)
        self.business_logic.generate_permanent_url = Mock(return_value=expected_permanent_url)

        # TODO: requires mocking of the S3 provider. implement later
        download_data = self.business_logic.get_dataset_artifact_download_data(dataset.version_id, artifact.id)
        expected_download_data = DatasetArtifactDownloadData(expected_file_size, expected_permanent_url)
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

    def test_get_ingest_manifest(self):
        published_version = self.initialize_published_collection()
        dataset = published_version.datasets[0]
        expected_manifest = IngestionManifest(
            anndata=[artifact.uri for artifact in dataset.artifacts if artifact.type == DatasetArtifactType.RAW_H5AD][
                0
            ],
            atac_fragment=[
                artifact.uri for artifact in dataset.artifacts if artifact.type == DatasetArtifactType.ATAC_FRAGMENT
            ][0],
        )
        manifest = self.business_logic.get_ingestion_manifest(dataset.version_id)
        self.assertEqual(expected_manifest, manifest)


class TestGetAllDatasets(BaseBusinessLogicTestCase):
    def test_get_all_private_datasets_ok(self):
        """
        Private datasets the user is authorized to view can be retrieved with
        `get_all_private_collection_versions_with_datasets`.
        """
        # test_user_1:
        # - private collection (2 datasets)
        # - public collection (2 datasets)
        # - published revision (2 datasets)
        # - unpublished revision, unchanged datasets (2 datasets)
        # - unpublished revision, new dataset, changed dataset, unchanged dataset (3 datasets)
        test_user_1 = "test_user_1"
        private_cv_1 = self.initialize_unpublished_collection(owner=test_user_1)
        self.initialize_published_collection(owner=test_user_1)
        self.initialize_collection_with_a_published_revision(owner=test_user_1)
        self.initialize_collection_with_an_unpublished_revision(owner=test_user_1)
        # Create unpublished revision with a replaced dataset and a new dataset.
        _, revision_1_updated = self.initialize_collection_with_an_unpublished_revision(owner=test_user_1)
        updated_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision_1_updated.version_id, "http://fake.url", None, revision_1_updated.datasets[0].version_id
        )
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision_1_updated.version_id, "http://fake.url", None, None
        )

        # test_user_2:
        # - private collection
        # - public collection
        # - published revision
        # - unpublished revision
        test_user_2 = "test_user_2"
        private_cv_2 = self.initialize_unpublished_collection(owner=test_user_2)
        self.initialize_published_collection(owner=test_user_2)
        self.initialize_collection_with_a_published_revision(owner=test_user_2)
        _, revision_2 = self.initialize_collection_with_an_unpublished_revision(owner=test_user_2)

        # Validate the expected datasets are returned.
        def _validate(actual: List[CollectionVersionWithDatasets], expected: List[CollectionVersionWithDatasets]):
            # Confirm the expected number of collection versions are returned.
            self.assertEqual(len(expected), len(actual))

            # Sort collection versions by ID, for comparison.
            actual.sort(key=lambda cv: cv.version_id.id)
            expected.sort(key=lambda cv: cv.version_id.id)

            # Confirm datasets length and content are correct for each collection version.
            for index, collection_version in enumerate(actual):
                datasets = collection_version.datasets
                expected_datasets = expected[index].datasets

                self.assertEqual(len(expected_datasets), len(datasets))

                # Check actual datasets match in expected.
                self.assertCountEqual(
                    [d.version_id for d in expected_datasets],
                    [d.version_id for d in datasets],
                )

        # Create the expected shape of revision_1_updated: datasets should only include the replacement dataset as
        # well as the new dataset.
        revision_1_updated_expected = deepcopy(revision_1_updated)
        revision_1_updated_expected.datasets = [
            self.database_provider.get_dataset_version(updated_dataset_version_id),
            self.database_provider.get_dataset_version(new_dataset_version_id),
        ]

        with self.subTest("With super user"):
            collection_versions = self.business_logic.get_private_collection_versions_with_datasets()
            # Super user should see private collections from both users, and the updated revision from test_user_1.
            expected = [private_cv_1, private_cv_2, revision_1_updated_expected]
            _validate(collection_versions, expected)

        with self.subTest("With owner"):
            collection_versions = self.business_logic.get_private_collection_versions_with_datasets(owner=test_user_1)
            # Owner should see their private collection, and their revision that has been udpated.
            expected = [private_cv_1, revision_1_updated_expected]
            _validate(collection_versions, expected)


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

    def test_update_dataset_status_validate_message_with_appending(self):
        """New messages are appended to the existing ones"""
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        dataset = unpublished_collection.datasets[0]
        error_message = "Validation error!"

        for _ in range(2):
            self.business_logic.update_dataset_version_status(
                dataset.version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID, error_message
            )
        version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
        validation_message = version_from_db.status.validation_message.split("\n")
        self.assertEqual([error_message] * 2, validation_message)

    def test_add_dataset_artifact_ok(self):
        """
        A dataset artifact can be added using `add_dataset_artifact`
        """
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)
        self.assertEqual(2, len(unpublished_collection.datasets))
        for dataset in unpublished_collection.datasets:
            self.assertEqual(dataset.artifacts, [])
            self.business_logic.add_dataset_artifact(
                dataset.version_id, DatasetArtifactType.H5AD, "http://fake.uri/artifact.h5ad"
            )

            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            self.assertEqual(1, len(version_from_db.artifacts))
            self.assertEqual(version_from_db.artifacts[0].type, DatasetArtifactType.H5AD.value)
            self.assertEqual(version_from_db.artifacts[0].uri, "http://fake.uri/artifact.h5ad")

    def test_update_dataset_artifact_ok(self):
        """
        A dataset artifact can be updated using `update_dataset_artifact`
        """
        published_collection = self.initialize_published_collection()
        for dataset in published_collection.datasets:
            cxg_artifact = [artifact for artifact in dataset.artifacts if artifact.type == "cxg"][0]
            self.assertEqual(cxg_artifact.uri, f"s3://cellxgene/{dataset.version_id}.cxg/")
            self.business_logic.update_dataset_artifact(cxg_artifact.id, "s3://cellxgene/new-name.cxg/")

            version_from_db = self.database_provider.get_dataset_version(dataset.version_id)
            updated_cxg_artifact = [
                artifact for artifact in version_from_db.artifacts if artifact.id == cxg_artifact.id
            ][0]
            self.assertEqual(updated_cxg_artifact.uri, "s3://cellxgene/new-name.cxg/")

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

    def test_update_dataset_artifact_metadata_missing_title_error(self):
        """
        Attempting to update a dataset without providing a title raises an exception.
        """

        # Create a private collection.
        unpublished_collection = self.initialize_unpublished_collection()

        # Create update.
        update = DatasetArtifactMetadataUpdate()

        # Confirm error is raised.
        dataset = unpublished_collection.datasets[0]
        with self.assertRaises(InvalidMetadataException):
            self.business_logic.update_dataset_artifact_metadata(
                unpublished_collection.version_id, dataset.version_id, update
            )

    def test_update_dataset_artifact_metadata_public_collection_error(self):
        """
        Attemping to update a dataset in a public collection raises an exception.
        """

        # Create public collection.
        published_collection = self.initialize_published_collection()

        # Create update.
        update = DatasetArtifactMetadataUpdate(title="new title")

        # Confirm error is raised.
        dataset = published_collection.datasets[0]
        with self.assertRaises(CollectionIsPublishedException):
            self.business_logic.update_dataset_artifact_metadata(
                published_collection.version_id, dataset.version_id, update
            )

    def test_update_dataset_artifact_metadata_invalid_status_error(self):
        """
        Attempting to update a dataset with a non-success processing status raises an exception.
        """

        # Create a private collection with a dataset that has not been fully ingested.
        unpublished_collection = self.initialize_unpublished_collection(complete_dataset_ingestion=False)

        # Create update.
        update = DatasetArtifactMetadataUpdate(title="new title")

        # Confirm error is raised.
        dataset = unpublished_collection.datasets[0]
        with self.assertRaises(DatasetInWrongStatusException):
            self.business_logic.update_dataset_artifact_metadata(
                unpublished_collection.version_id, dataset.version_id, update
            )

    def test_update_dataset_artifact_metadata_ok(self):
        """
        The dataset artifact metadata can be updated using `update_dataset_artifact_metadata`
        """

        # Create a private collection.
        unpublished_collection = self.initialize_unpublished_collection()

        # Create update.
        update = DatasetArtifactMetadataUpdate(title="new title")

        # Mock trigger artifact update.
        self.business_logic.trigger_dataset_artifact_update = Mock()

        # Attempt update.
        dataset = unpublished_collection.datasets[0]
        self.business_logic.update_dataset_artifact_metadata(
            unpublished_collection.version_id, dataset.version_id, update
        )

        # Confirm trigger was called.
        self.business_logic.trigger_dataset_artifact_update.assert_called_once()

    def test_update_dataset_artifact_metadata_sanitized_ok(self):
        """
        The dataset artifact metadata can be updated with sanitized values.
        """

        # Create a private collection.
        unpublished_collection = self.initialize_unpublished_collection()

        # Create update.
        title = "new title"
        update = DatasetArtifactMetadataUpdate(title=f" {title} ")

        # Mock trigger artifact update.
        trigger_mock = self.business_logic.trigger_dataset_artifact_update = Mock()

        # Attempt update.
        dataset = unpublished_collection.datasets[0]
        self.business_logic.update_dataset_artifact_metadata(
            unpublished_collection.version_id, dataset.version_id, update
        )

        # Confirm update was sanitized.
        args, _ = trigger_mock.call_args
        assert args[1].title == title

    def test_update_dataset_artifact_metadata_invalid_title_error(self):
        """
        Attempting to update a dataset with an invalid title raises an exception.
        """

        # Create a private collection.
        unpublished_collection = self.initialize_unpublished_collection()

        # Create update.
        update = DatasetArtifactMetadataUpdate(title="new\x00title")

        # Confirm error is raised.
        dataset = unpublished_collection.datasets[0]
        with self.assertRaises(InvalidMetadataException):
            self.business_logic.update_dataset_artifact_metadata(
                unpublished_collection.version_id, dataset.version_id, update
            )


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

    def test_create_collection_version_multiple_auto_versions_false_fails(self):
        """
        Given a canonical collection with a published version and an unpublished version with is_auto_version == False,
        creating an unpublished collection version where is_auto_verison is False will fail.
        """
        published_collection = self.initialize_published_collection()
        self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=False)

        with self.assertRaises(CollectionVersionException):
            self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=False)

    def test_create_collection_version_multiple_auto_versions_fails(self):
        """
        Given a canonical collection with a published version and an unpublished version with is_auto_version == True,
        creating another unpublished collection version where is_auto_verison is False will fail.
        """
        published_collection = self.initialize_published_collection()
        self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=True)

        with self.assertRaises(CollectionVersionException):
            self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=False)

    def test_create_collection_version_multiple_auto_versions_true_fails(self):
        """
        Given a canonical collection with a published version and an unpublished version where is_auto_version is True,
        creating another unpublished collection version where is_auto_version is True will fail.
        """
        published_collection = self.initialize_published_collection()
        self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=True)

        with self.assertRaises(CollectionVersionException):
            self.business_logic.create_collection_version(published_collection.collection_id, is_auto_version=True)

    def test_create_collection_version_is_auto_version_succeeds(self):
        """
        Given a canonical collection with a published version and no unpublished versions,
        creating an unpublished collection version where is_auto_verison is True will succeed.
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=True
        )

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

        self.assertTrue(new_version.is_auto_version)

    def test_create_collection_version_is_auto_version_succeeds__with_revision(self):
        """
        Given a canonical collection with a published version and an unpublished version where is_auto_version is False,
        creating another unpublished collection version where is_auto_version is True will succeed.
        """
        published_collection = self.initialize_published_collection()
        existing_revision = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=False
        )
        auto_revision = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=True
        )

        # The new version has a different version_id
        self.assertNotEqual(published_collection.version_id, auto_revision.version_id)
        self.assertNotEqual(existing_revision.version_id, auto_revision.version_id)

        # Both revisions are linked to the same canonical collection
        self.assertEqual(published_collection.collection_id, existing_revision.collection_id)
        self.assertEqual(published_collection.collection_id, auto_revision.collection_id)

        # The new version has the same collection metadata, and datasets
        self.assertEqual(published_collection.metadata, auto_revision.metadata)
        self.assertEqual(published_collection.datasets, auto_revision.datasets)

        # The new version should not be published (aka Private)
        self.assertIsNone(auto_revision.published_at)

        # The new version has is_auto_version == True, while existing is False
        self.assertFalse(existing_revision.is_auto_version)
        self.assertTrue(auto_revision.is_auto_version)

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

        self.business_logic.delete_collection_version(new_version)

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
        A Collection can be tombstoned
        """
        published_collection = self.initialize_published_collection()

        public_dataset_asset_s3_uris = [
            f"s3://datasets/{dv.version_id}.{ext}" for ext in ("rds", "h5ad") for dv in published_collection.datasets
        ]

        # Verify public Dataset asset files are in place in s3 store
        assert all(self.s3_provider.uri_exists(uri) for uri in public_dataset_asset_s3_uris)

        self.business_logic.tombstone_collection(published_collection.collection_id)

        # The collection version canonical collection has tombstoned attribute marked as True
        collection = self.business_logic.get_canonical_collection(published_collection.collection_id)
        assert collection.tombstoned is True
        # Verify public Dataset asset files are gone
        assert all(not self.s3_provider.uri_exists(uri) for uri in public_dataset_asset_s3_uris)

    def test_resurrect_collection_ok(self):
        """
        A tombstoned Collection can be resurrected
        """
        published_collection = self.initialize_published_collection()
        public_dataset_asset_s3_uris = [
            f"s3://datasets/{dv.version_id}.{ext}" for ext in ("rds", "h5ad") for dv in published_collection.datasets
        ]

        self.business_logic.tombstone_collection(published_collection.collection_id)
        self.business_logic.resurrect_collection(published_collection.collection_id)

        # The collection is no longer tombstoned
        canonical_collection = self.business_logic.get_canonical_collection(published_collection.collection_id)
        assert canonical_collection.tombstoned is False
        # Verify public Dataset asset files are restored
        assert all(self.s3_provider.uri_exists(uri) for uri in public_dataset_asset_s3_uris)

    def test_resurrect_collection_with_tombstoned_dataset_ok(self):
        """
        Individually-tombstoned Datasets should remain tombstoned after resurrection of their parent Collection.
        """
        published_collection = self.initialize_published_collection()
        dataset_to_keep, dataset_to_tombstone = published_collection.datasets
        public_dataset_asset_s3_uris_kept = [
            f"s3://datasets/{dataset_to_keep.version_id}.{ext}" for ext in ("rds", "h5ad")
        ]
        public_dataset_asset_s3_uris_tombstoned = [
            f"s3://datasets/{dataset_to_tombstone.version_id}.{ext}" for ext in ("rds", "h5ad")
        ]

        # Tombstone the Dataset
        revision = self.business_logic.create_collection_version(published_collection.collection_id)
        self.business_logic.remove_dataset_version(revision.version_id, dataset_to_tombstone.version_id, True)
        self.business_logic.publish_collection_version(revision.version_id)

        self.business_logic.tombstone_collection(published_collection.collection_id)
        self.business_logic.resurrect_collection(published_collection.collection_id)

        # The Collection is no longer tombstoned
        canonical_collection = self.business_logic.get_canonical_collection(published_collection.collection_id)
        assert canonical_collection.tombstoned is False

        # Dataset that was kept should not be tombstoned
        dataset_kept = self.business_logic.get_dataset_version_from_canonical(
            dataset_to_keep.dataset_id, get_tombstoned=True
        )
        assert dataset_kept.canonical_dataset.tombstoned is False
        # Dataset that was tombstoned should still be tombstoned
        dataset_is_tombstoned_exception_raised = False
        try:
            self.business_logic.get_dataset_version_from_canonical(dataset_to_tombstone.dataset_id, get_tombstoned=True)
        except DatasetIsTombstonedException:
            dataset_is_tombstoned_exception_raised = True
        assert dataset_is_tombstoned_exception_raised

        # Public-access Dataset asset files for kept Dataset are restored
        assert all(self.s3_provider.uri_exists(uri) for uri in public_dataset_asset_s3_uris_kept)
        # Public-access Dataset asset files for tombstoned Dataset are not restored
        assert all(not self.s3_provider.uri_exists(uri) for uri in public_dataset_asset_s3_uris_tombstoned)

    def test_publish_version_fails_on_published_collection(self):
        """
        `publish_collection_version` should fail if called on a collection version that is already published.
        """
        published_collection = self.initialize_published_collection()
        self.assertIsNotNone(published_collection.published_at)

        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(published_collection.version_id)

    def test_publish_version_fails_on_mismatching_dataset_schema_versions(self):
        """
        `publish_collection_version` should fail if called on a collection version with datasets that have mismatching
        schema versions
        """
        cv = self.initialize_unpublished_collection()
        metadata = deepcopy(self.sample_dataset_metadata)
        metadata.schema_version = "3.1.0"
        self.database_provider.set_dataset_metadata(cv.datasets[0].version_id, metadata)

        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(cv.version_id)

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
        self.assertEqual(published_version.data_submission_policy_version, DATA_SUBMISSION_POLICY_VERSION)
        self.assertEqual(published_version.is_auto_version, False)

        # The published and unpublished collection have the same collection_id and version_id
        self.assertEqual(published_version.collection_id, unpublished_collection.collection_id)
        self.assertEqual(published_version.version_id, unpublished_collection.version_id)

        # get_collection retrieves the correct version
        version = self.business_logic.get_published_collection_version(unpublished_collection.collection_id)
        if version:  # pylance
            self.assertEqual(version.version_id, published_version.version_id)

    def test_collection_url(self):
        self.assertEqual(
            self.business_logic.get_collection_url(collection_id="test-collection-id"),
            "https://collections_domain/collections/test-collection-id",
        )

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

        self.business_logic.remove_dataset_version(
            new_version.version_id, dataset_version_to_remove.version_id, delete_published=True
        )

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
        self.business_logic.set_dataset_metadata(added_dataset_version_id, self.sample_dataset_metadata)
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
        self.business_logic.set_dataset_metadata(replaced_dataset_version_id, self.sample_dataset_metadata)
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

    def test_publish_version_sets_is_auto_version_False(self):
        """
        When publishing a collection version, is_auto_version should be set to False
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=True
        )
        self.assertTrue(new_version.is_auto_version)
        self.business_logic.publish_collection_version(new_version.version_id)

        published_new_version = self.business_logic.get_collection_version(new_version.version_id)
        self.assertFalse(published_new_version.is_auto_version)

    def test_publish_version_cleanup_with_auto_version(self):
        """
        When publishing an auto_version, do NOT cleanup dataset versions from an open revision
        """
        published_collection = self.initialize_published_collection()
        collection_open_revision = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=False
        )
        dataset_to_replace = collection_open_revision.datasets[0]
        new_dataset_version_in_open_revision = self.database_provider.replace_dataset_in_collection_version(
            collection_open_revision.version_id, dataset_to_replace.version_id
        )
        auto_version = self.business_logic.create_collection_version(
            published_collection.collection_id, is_auto_version=True
        )
        self.business_logic.publish_collection_version(auto_version.version_id)

        # open revision's new dataset and its artifacts should still be persisted
        self.assertIsNotNone(
            self.database_provider.get_dataset_version(new_dataset_version_in_open_revision.version_id)
        )
        [self.assertTrue(self.s3_provider.uri_exists(a.uri)) for a in dataset_to_replace.artifacts]

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
        self.business_logic.set_dataset_metadata(added_dataset_version_id, self.sample_dataset_metadata)
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
        self.business_logic.set_dataset_metadata(replaced_dataset_version_id, self.sample_dataset_metadata)
        self.business_logic.publish_collection_version(new_version.version_id)

        version_from_db = self.database_provider.get_collection_version(new_version.version_id)
        dataset_0 = self.database_provider.get_dataset_version(version_from_db.datasets[0])
        dataset_1 = self.database_provider.get_dataset_version(version_from_db.datasets[1])

        self.assertNotEqual(dataset_0.version_id, dataset_id_to_replace)
        self.assertEqual(dataset_1.version_id, dataset_id_to_keep)

        self.assertEqual(dataset_0.canonical_dataset.published_at, published_collection.published_at)
        self.assertEqual(dataset_1.canonical_dataset.published_at, published_collection.published_at)

    def test_publish_new_collection_no_doi_crossref_not_checked_ok(self):
        """
        When publishing a new collection without a DOI, no Crossref check should be made.
        """

        # Create an unpublished collection without a DOI.
        collection = self.initialize_unpublished_collection()

        # Mock the Crossref check.
        self.crossref_provider.fetch_metadata = Mock()

        # Publish the collection.
        self.business_logic.publish_collection_version(collection.version_id)

        # Confirm Crossref was not called.
        self.crossref_provider.fetch_metadata.assert_not_called()

    def test_publish_private_collection_crossref_checked_ok(self):
        """
        When publishing a new collection with a DOI, a Crossref check should be made.
        """

        # Mock the Crossref check.
        doi_curie = "test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, doi_curie, 17169328.664)
        )

        # Create a private collection with a DOI.
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        collection_version = self.initialize_unpublished_collection(links=links)

        # Publish the private collection.
        self.business_logic.publish_collection_version(collection_version.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

    def test_publish_private_collection_same_doi_publisher_metadata_not_updated_ok(self):
        """
        When publishing a private collection - and there is no change in DOI - publisher metadata should not
        be updated if the deposited date is before the created at date of the corresponding collection.
        """

        # Get the time before the collection is created.
        before_created_at = datetime.utcnow().timestamp()

        # Create a private collection with a DOI.
        doi_curie = "test/doi"
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        collection_version = self.initialize_unpublished_collection(links=links)

        # Mock the Crossref update, and set the deposited date to be before the created at date of the collection.
        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, doi_curie, before_created_at)
        )

        # Mock the publisher metadata update.
        self.database_provider.save_collection_publisher_metadata = Mock()

        # Publish the private collection.
        self.business_logic.publish_collection_version(collection_version.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

        # Confirm the metadata was not updated.
        self.database_provider.save_collection_publisher_metadata.assert_not_called()

    def test_publish_private_collection_same_doi_publisher_metadata_updated_ok(self):
        """
        When publishing a private collection - and there is no change in DOI - publisher metadata should be
        updated if the deposited date is after the created at date of the corresponding collection.
        """

        # Create a private collection.
        doi_curie = "test/doi"
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        collection_version = self.initialize_unpublished_collection(links=links)

        # Get the time after the collection has been created.
        after_created_at = datetime.utcnow().timestamp()

        # Mock the Crossref update, and set the deposited date to be after the created at date of the collection.
        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, doi_curie, after_created_at)
        )

        # Mock the publisher metadata update.
        self.database_provider.save_collection_publisher_metadata = Mock()

        # Publish the private collection.
        self.business_logic.publish_collection_version(collection_version.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

        # Confirm the publisher metadata was updated.
        self.database_provider.save_collection_publisher_metadata.assert_called_once()

    def test_publish_revision_no_doi_crossref_not_checked_ok(self):
        """
        Crossref should not be called to check for updates to the publisher metadata when publishing
        a revision that has no DOI.
        """

        # Create a revision without a DOI.
        _, revision = self.initialize_collection_with_an_unpublished_revision()

        # Mock the Crossref check.
        self.crossref_provider.fetch_metadata = Mock()

        # Publish the revision.
        self.business_logic.publish_collection_version(revision.version_id)

        # Confirm Crossref was not called.
        self.crossref_provider.fetch_metadata.assert_not_called()

    def test_publish_revision_crossref_checked_ok(self):
        """
        Crossref should be called to check for updates to the publisher metadata when publishing
        a revision that has a DOI.
        """

        # Mock the Crossref check.
        doi_curie = "test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, doi_curie, 17169328.664)
        )

        # Create a revision with a DOI.
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        _, revision = self.initialize_collection_with_an_unpublished_revision(links=links)

        # Reset the Crossref check so we can check it on publish.
        self.crossref_provider.fetch_metadata.reset_mock()

        # Publish the revision.
        self.business_logic.publish_collection_version(revision.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

    def test_publish_revision_changed_doi_collection_updated_ok(self):
        """
        When publishing a revision - and there is a change in DOI - collection and artifacts should be updated.
        """

        # Mock the Crossref check.
        doi_curie = "test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, doi_curie, 17169328.664)
        )

        # Create a revision with a DOI.
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        _, revision = self.initialize_collection_with_an_unpublished_revision(links=links)

        # Mock a DOI change from Crossref.
        new_doi_curie = "new/test/doi"
        self.crossref_provider.fetch_metadata.reset_mock()
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, new_doi_curie, 17169328.664)
        )

        # Mock the call to update the collection and artifacts.
        self.business_logic.update_collection_version = Mock()

        # Mock the publisher metadata update.
        self.database_provider.save_collection_publisher_metadata = Mock()

        # Publish the revision.
        with self.assertRaises(CollectionPublishException):
            self.business_logic.publish_collection_version(revision.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

        # Confirm collection (and artifacts) update was called.
        self.business_logic.update_collection_version.assert_called_once()

        # Confirm the metadata was not updated.
        self.database_provider.save_collection_publisher_metadata.assert_not_called()

    def test_publish_revision_same_doi_publisher_metadata_not_updated_ok(self):
        """
        When publishing a revision - and there is no change in DOI - publisher metadata should not be
        updated if the deposited date is before the revised at or originally published at date of the
        corresponding collection.
        """

        # Get the time before the collection is published.
        before_published_at = datetime.utcnow().timestamp()

        # Create a revision with a DOI.
        doi_curie = "test/doi"
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        _, revision = self.initialize_collection_with_an_unpublished_revision(links=links)

        # Mock the Crossref update, and set the deposited date to be before the published date of the collection.
        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, doi_curie, before_published_at)
        )

        # Mock the publisher metadata update.
        self.database_provider.save_collection_publisher_metadata = Mock()

        # Publish the revision.
        self.business_logic.publish_collection_version(revision.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

        # Confirm the metadata was not updated.
        self.database_provider.save_collection_publisher_metadata.assert_not_called()

    def test_publish_revision_same_doi_publisher_metadata_updated_ok(self):
        """
        When publishing a revision - and there is no change in DOI - publisher metadata should be updated
        if the deposited date is after the revised at or published at date of the corresponding collection.
        """

        # Create a published collection.
        doi_curie = "test/doi"
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        collection = self.initialize_published_collection(links=links)

        # Get the time after the collection has been published.
        after_published_at = datetime.utcnow().timestamp()

        # Create a revision.
        revision = self.business_logic.create_collection_version(collection.collection_id)

        # Mock the Crossref update, and set the deposited date to be after the published date of the collection.
        expected_publisher_metadata = {"authors": ["Test Author"]}
        self.crossref_provider.fetch_metadata = Mock(
            return_value=(expected_publisher_metadata, doi_curie, after_published_at)
        )

        # Mock the publisher metadata update.
        self.database_provider.save_collection_publisher_metadata = Mock()

        # Publish the revision.
        self.business_logic.publish_collection_version(revision.version_id)

        # Confirm Crossref was called.
        self.crossref_provider.fetch_metadata.assert_called_once()

        # Confirm the publisher metadata was updated.
        self.database_provider.save_collection_publisher_metadata.assert_called_once()

    def test_update_crossref_metadata_same_doi_returns_none_ok(self):
        """
        The check for changes in Crossref returns None if there is no change between the existing DOI and the
        DOI returned from Crossref.
        """

        # Mock the Crossref check.
        doi_curie = "test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, doi_curie, 17169328.664)
        )

        # Create a revision with a DOI.
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        _, revision = self.initialize_collection_with_an_unpublished_revision(links=links)

        # Check Crossref for updates.
        doi_update = self.business_logic._update_crossref_metadata(revision, datetime.utcnow())

        # Confirm the returned DOI update is None.
        self.assertIsNone(doi_update)

    def test_update_crossref_metadata_changed_doi_returns_doi_update_ok(self):
        """
        The check for changes in Crossref returns the existing DOI and the returned Crossref DOI if there is a
        change between the two.
        """

        # Mock the Crossref check.
        doi_curie = "test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, doi_curie, 17169328.664)
        )

        # Create a revision with a DOI.
        links = [Link("test doi", "DOI", f"https://doi.org/{doi_curie}")]
        _, revision = self.initialize_collection_with_an_unpublished_revision(links=links)

        # Mock a DOI change from Crossref.
        new_doi_curie = "new/test/doi"
        self.crossref_provider.fetch_metadata = Mock(
            return_value=({"authors": ["Test Author"]}, new_doi_curie, 17169328.664)
        )

        # Check Crossref for updates.
        doi_update = self.business_logic._update_crossref_metadata(revision, datetime.utcnow())

        # Confirm the DOI update is not None, and the existing and new DOI are returned.
        self.assertIsNotNone(doi_update)
        self.assertEqual(doi_update[0], f"https://doi.org/{doi_curie}")
        self.assertEqual(doi_update[1], f"https://doi.org/{new_doi_curie}")


class TestCollectionUtilities(BaseBusinessLogicTestCase):
    def test__delete_datasets_from_public_access_bucket(self):
        """
        Test Dataset deletes when tombstoning a public Collection with published and updated Datasets in an outstanding
        Revision
        """
        published_collection = self.initialize_published_collection()
        new_version = self.business_logic.create_collection_version(published_collection.collection_id)

        # Update the first dataset
        dataset_id_to_replace = published_collection.datasets[0].version_id

        replaced_dataset_version_id, _ = self.business_logic.ingest_dataset(
            new_version.version_id, "http://fake.url", None, dataset_id_to_replace
        )

        self.complete_dataset_processing_with_success(replaced_dataset_version_id)

        dataset_versions = published_collection.datasets + [
            self.business_logic.get_dataset_version(replaced_dataset_version_id)
        ]
        expected_delete_keys = set()
        fake_public_bucket = "datasets"
        for d_v in dataset_versions:
            for file_type in ("h5ad", "rds"):
                key = f"{d_v.version_id}.{file_type}"
                self.s3_provider.upload_file(None, fake_public_bucket, key, None)  # Populate s3 mock with assets
                self.assertTrue(self.s3_provider.uri_exists(f"s3://{fake_public_bucket}/{key}"))
                expected_delete_keys.add(f"{d_v.version_id}.{file_type}")
            expected_delete_keys.update(self.business_logic.get_atac_fragment_uris_from_dataset_version(d_v))
        self.assertTrue(len(expected_delete_keys) > 0)
        self.assertTrue(all(self.s3_provider.file_exists(fake_public_bucket, key) for key in expected_delete_keys))
        actual_delete_keys = set(
            self.business_logic.delete_all_dataset_versions_from_public_bucket_for_collection(
                published_collection.collection_id
            )
        )
        self.assertEqual(expected_delete_keys, actual_delete_keys)

    def test__restore_previous_collection_version(self):
        """
        Test restoring a previous version of a collection
        """
        original_collection_version = self.initialize_published_collection()
        collection_id = original_collection_version.collection_id
        new_version = self.business_logic.create_collection_version(collection_id)
        self.business_logic.publish_collection_version(new_version.version_id)

        # Restore the previous version
        self.business_logic.restore_previous_collection_version(collection_id)

        # Fetch the collection
        restored_collection = self.business_logic.get_canonical_collection(collection_id)

        # Ensure the collection is pointing back at the original collection version
        self.assertEqual(original_collection_version.version_id, restored_collection.version_id)

        # Ensure replaced CollectionVersion is deleted from DB
        self.assertIsNone(self.database_provider.get_collection_version(new_version.version_id))

    def test__restore_previous_collection_version__no_previous_versions(self):
        """
        Test restoring a previous version of a collection fails when there are no previous versions
        """
        original_collection_version = self.initialize_published_collection()
        collection_id = original_collection_version.collection_id

        with self.assertRaises(NoPreviousCollectionVersionException):
            self.business_logic.restore_previous_collection_version(collection_id)


class TestGetEntitiesBySchema(BaseBusinessLogicTestCase):
    def test_get_latest_published_collection_versions_by_schema(self):
        """
        Test fetching index of most recently published collection versions with datasets matching a given
        CxG schema version
        """
        # set-up: 1 private, 3 public datasets
        # public just 3.1.0
        cv = self.initialize_unpublished_collection(num_datasets=1)
        metadata = deepcopy(self.sample_dataset_metadata)
        metadata.schema_version = "3.1.0"
        self.database_provider.set_dataset_metadata(cv.datasets[0].version_id, metadata)
        self.business_logic.publish_collection_version(cv.version_id)

        # public with 4.0.0, and an unpublished revision with 4.1.0 (should not appear)
        cv = self.initialize_unpublished_collection(num_datasets=1)
        metadata = deepcopy(self.sample_dataset_metadata)
        metadata.schema_version = "4.0.0"
        self.database_provider.set_dataset_metadata(cv.datasets[0].version_id, metadata)
        self.business_logic.publish_collection_version(cv.version_id)

        revision = self.business_logic.create_collection_version(cv.collection_id)
        metadata.schema_version = "4.1.0"
        self.database_provider.set_dataset_metadata(revision.datasets[0].version_id, metadata)

        # public with 3.0.0 and 3.1.0
        _, cv = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        metadata = deepcopy(self.sample_dataset_metadata)
        metadata.schema_version = "3.1.0"
        self.database_provider.set_dataset_metadata(cv.datasets[0].version_id, metadata)
        self.business_logic.publish_collection_version(cv.version_id)

        # private with just 4.1.0, should not appear
        cv = self.initialize_unpublished_collection(num_datasets=1)
        metadata = deepcopy(self.sample_dataset_metadata)
        metadata.schema_version = "4.1.0"
        self.database_provider.set_dataset_metadata(cv.datasets[0].version_id, metadata)

        with self.subTest("Without schema_version Parameter"):
            collection_versions = self.business_logic.get_all_mapped_collection_versions_with_datasets()
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["3.1.0", "3.1.0", "4.0.0"])
        with self.subTest("Querying against major schema version 3"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("3._._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["3.1.0", "3.1.0"])
        with self.subTest("Querying against major schema version 4"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("4._._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["4.0.0"])
        with self.subTest("Querying against minor schema version 3.0"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("3.0._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["3.0.0"])
        with self.subTest("Querying against minor schema version 3.1"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("3.1._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["3.1.0", "3.1.0"])
        with self.subTest("Querying against patch schema version 4.0.0"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("4.0.0")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], ["4.0.0"])
        with self.subTest("Querying against patch schema version 4.1.0"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("4.1.0")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], [])
        with self.subTest("Querying against minor schema version 4.1"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("4.1._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], [])
        with self.subTest("Querying against major schema version 5"):
            collection_versions = self.business_logic.get_latest_published_collection_versions_by_schema("5._._")
            self.assertCountEqual([cv.schema_version for cv in collection_versions], [])


class TestSetPreviousDatasetVersion(BaseBusinessLogicTestCase):
    def test_with_published_collection(self):
        collection = self.initialize_published_collection(num_datasets=1)
        dataset = collection.datasets[0]
        with self.assertRaises(CollectionIsPublishedException):
            self.business_logic.restore_previous_dataset_version(collection.version_id, dataset.dataset_id)

    def test_private_with_no_previous_version(self):
        collection = self.initialize_unpublished_collection(num_datasets=1)
        dataset = collection.datasets[0]
        with self.assertRaises(NoPreviousDatasetVersionException):
            self.business_logic.restore_previous_dataset_version(collection.version_id, dataset.dataset_id)

    def test_revision_with_no_previous_version(self):
        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        dataset = revision.datasets[0]
        with self.assertRaises(NoPreviousDatasetVersionException):
            self.business_logic.restore_previous_dataset_version(revision.version_id, dataset.dataset_id)

    def test_private_restoring_previously_replaced_dataset(self):
        collection = self.initialize_unpublished_collection(num_datasets=1)
        previous_dataset = collection.datasets[0]
        dataset_id = previous_dataset.dataset_id
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection.version_id, "http://fake.url", None, previous_dataset.version_id
        )

        # Verify the new dataset is in the collection and the previous dataset is not
        datasets = [d.version_id for d in self.business_logic.get_collection_version(collection.version_id).datasets]
        assert new_dataset_version_id in datasets
        assert previous_dataset.version_id not in datasets

        # Restore the previous dataset
        self.business_logic.restore_previous_dataset_version(collection.version_id, dataset_id)

        # Verify the new dataset is no longer in the collection and the previous one is restored
        datasets = [d.version_id for d in self.business_logic.get_collection_version(collection.version_id).datasets]
        assert new_dataset_version_id not in datasets
        assert previous_dataset.version_id in datasets

    def test_revision_restoring_previously_published_dataset(self):
        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        previous_dataset = revision.datasets[0]
        dataset_id = previous_dataset.dataset_id
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id, "http://fake.url", None, previous_dataset.version_id
        )

        # Verify the new dataset is in the collection and the previous dataset is not
        datasets = [d.version_id for d in self.business_logic.get_collection_version(revision.version_id).datasets]
        assert new_dataset_version_id in datasets
        assert previous_dataset.version_id not in datasets

        # Restore the previous dataset
        self.business_logic.restore_previous_dataset_version(revision.version_id, dataset_id)

        # Verify the new dataset is no longer in the collection and the previous one is restored
        datasets = [d.version_id for d in self.business_logic.get_collection_version(revision.version_id).datasets]
        assert new_dataset_version_id not in datasets
        assert previous_dataset.version_id in datasets


class TestConcurrentUpdates(BaseBusinessLogicTestCase):
    """
    Test that concurrent updates to a "shared" field in DatasetVersion are properly supported.
    This is simulated by calling a method that updates such field concurrently, using a thread pool.
    Note: this test should always pass, but we can't guarantee that it actually prevents a race condition.
    Different hosts might have different behaviors and might not yield a race condition in the first place
    (for instance, if the runtime can only give 1 worker to ThreadPoolExecutor).
    """

    def test_concurrency(self):
        collection = self.initialize_published_collection(num_datasets=1)
        dataset = collection.datasets[0]

        def add_artifact():
            self.database_provider.create_dataset_artifact(dataset.version_id, DatasetArtifactType.H5AD, "fake_uri")

        self.assertEqual(len(dataset.artifacts), 6)

        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=4) as e:
            for _ in range(10):
                e.submit(add_artifact)

        dv = self.business_logic.get_dataset_version(dataset.version_id)
        self.assertIsNotNone(dv)
        if dv is not None:
            self.assertEqual(len(dv.artifacts), 16)


class TestDatasetArtifactMetadataUpdates(BaseBusinessLogicTestCase):
    """
    Test methods related to updating metadata fields on dataset artifacts
    """

    def test_generate_dataset_citation(self):
        collection = self.initialize_unpublished_collection(num_datasets=1)
        dataset_version_id = collection.datasets[0].version_id
        doi = "https://doi.org/12.2345/science.abc1234"
        expected = (
            f"Publication: {doi} Dataset Version: {self.mock_config.dataset_assets_base_url}/{dataset_version_id}.h5ad "
            "curated and distributed by CZ CELLxGENE Discover in Collection: "
            f"{self.mock_config.collections_base_url}/collections/{collection.collection_id}"
        )

        self.assertEqual(
            expected,
            self.business_logic.generate_dataset_citation(collection.collection_id, dataset_version_id, doi),
        )

    def test_generate_dataset_citation__no_doi(self):
        collection = self.initialize_unpublished_collection(num_datasets=1)
        dataset_version_id = collection.datasets[0].version_id
        expected = (
            f"Dataset Version: {self.mock_config.dataset_assets_base_url}/{dataset_version_id}.h5ad "
            "curated and distributed by CZ CELLxGENE Discover in Collection: "
            f"{self.mock_config.collections_base_url}/collections/{collection.collection_id}"
        )

        self.assertEqual(
            expected, self.business_logic.generate_dataset_citation(collection.collection_id, dataset_version_id)
        )

    def test_trigger_dataset_artifact_update(self):
        metadata_update = DatasetArtifactMetadataUpdate(schema_version="4.0.0")

        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        current_dataset_version_id = revision.datasets[0].version_id
        self.batch_job_provider.start_metadata_update_batch_job = Mock()
        self.business_logic.trigger_dataset_artifact_update(revision, metadata_update, current_dataset_version_id)

        self.batch_job_provider.start_metadata_update_batch_job.assert_called_once()
        # Confirm citation is updated if new dataset version is generated as part of trigger_dataset_artifact_update
        self.assertIsNotNone(metadata_update.citation)

    def test_trigger_dataset_artifact_update__with_new_dataset_version_id(self):
        metadata_update = DatasetArtifactMetadataUpdate(schema_version="4.0.0")

        _, revision = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        current_dataset_version_id = revision.datasets[0].version_id
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            revision.version_id,
            "http://fake.url",
            None,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )

        self.batch_job_provider.start_metadata_update_batch_job = Mock()
        self.business_logic.trigger_dataset_artifact_update(
            revision, metadata_update, current_dataset_version_id, new_dataset_version_id
        )
        self.batch_job_provider.start_metadata_update_batch_job.assert_called_once_with(
            current_dataset_version_id, new_dataset_version_id, metadata_update
        )

    def test_trigger_dataset_artifact_update__error_if_published(self):
        metadata_update = DatasetArtifactMetadataUpdate(schema_version="4.0.0")
        collection, _ = self.initialize_collection_with_an_unpublished_revision(num_datasets=1)
        current_dataset_version_id = collection.datasets[0].version_id
        with self.assertRaises(CollectionIsPublishedException):
            self.business_logic.trigger_dataset_artifact_update(collection, metadata_update, current_dataset_version_id)


if __name__ == "__main__":
    unittest.main()
