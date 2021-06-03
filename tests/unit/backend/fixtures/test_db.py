from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    CollectionLinkType,
    DatasetArtifactType,
    DatasetArtifactFileType,
    DbCollection,
    DbCollectionLink,
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    UploadStatus,
    ValidationStatus,
    ConversionStatus,
    ProcessingStatus,
    DbGeneset,
)
from backend.corpora.common.utils.db_session import DBSessionMaker
from backend.corpora.common.utils.s3_buckets import cxg_bucket
from backend.scripts.create_db import create_db
from tests.unit.backend.fixtures import config


class TestDatabaseManager:
    is_initialized = False

    @classmethod
    def initialize_db(cls):
        if cls.is_initialized:
            return
        testdb = TestDatabase()
        testdb.create_db()
        testdb.populate_test_data()
        cls.is_initialized = True


class TestDatabase:
    fake_s3_file = f"s3://{config.CORPORA_TEST_CONFIG['bucket_name']}/test_s3_uri.h5ad"
    real_s3_file = "s3://corpora-data-dev/fake-h5ad-file.h5ad"

    def __init__(self, real_data=False):
        self.real_data = real_data

    def create_db(self):
        create_db()

    def populate_test_data(self):
        self.session = DBSessionMaker().session()
        self._populate_test_data()
        self.session.close()
        del self.session

    def _populate_test_data(self):
        self._create_test_collections()
        self._create_test_collection_links()
        self._create_test_datasets()
        self._create_test_dataset_artifacts()
        self._create_test_dataset_processing_status()
        self._create_test_geneset()

    def _create_test_collections(self):
        collection = DbCollection(
            id="test_collection_id",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
        )
        self.session.add(collection)
        collection = DbCollection(
            id="test_collection_id",
            visibility=CollectionVisibility.PRIVATE.name,
            owner="test_user_id",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
        )
        self.session.add(collection)
        collection = DbCollection(
            id="test_collection_id_public",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            name="test_collection_id_public",
            description="test_description",
            data_submission_policy_version="0",
        )
        self.session.add(collection)
        collection = DbCollection(
            id="test_collection_id_not_owner",
            visibility=CollectionVisibility.PRIVATE.name,
            owner="Someone_else",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
        )
        self.session.add(collection)
        self.session.commit()

    def _create_test_geneset(self):
        geneset = DbGeneset(
            id="test_geneset",
            name="test_geneset",
            description="this is a geneset",
            collection_id="test_collection_id",
            collection_visibility=CollectionVisibility.PUBLIC.name,
        )

        self.session.add(geneset)
        dataset = self.session.query(DbDataset).get("test_dataset_id")
        geneset = DbGeneset(
            id="test_geneset_with_dataset",
            name="test_geneset_with_dataset",
            description="this is a geneset with a dataset",
            collection_id="test_collection_id",
            collection_visibility=CollectionVisibility.PUBLIC.name,
            datasets=[dataset],
        )
        self.session.add(geneset)
        self.session.commit()

    def _create_test_collection_links(self):
        for link_type in CollectionLinkType:
            self.session.add(
                DbCollectionLink(
                    id=f"test_collection_{link_type.value}_link_id",
                    collection_id="test_collection_id",
                    collection_visibility=CollectionVisibility.PUBLIC.name,
                    link_name=f"test_{link_type.value}_link_name",
                    link_url=f"http://test_{link_type.value}_url.place",
                    link_type=link_type.name,
                )
            )
            self.session.add(
                DbCollectionLink(
                    id=f"test_collection_no_name_{link_type.value}_link_id",
                    collection_id="test_collection_id",
                    collection_visibility=CollectionVisibility.PUBLIC.name,
                    link_url=f"http://test_no_link_name_{link_type.value}_url.place",
                    link_type=link_type.name,
                )
            )
        self.session.commit()

    def _create_test_datasets(self):
        test_dataset_id = "test_dataset_id"
        dataset = DbDataset(
            id=test_dataset_id,
            revision=0,
            name="test_dataset_name",
            organism={"ontology_term_id": "test_obo", "label": "test_organism"},
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=["test_sex", "test_sex2"],
            ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            collection_id="test_collection_id",
            explorer_url="test_url",
            explorer_s3_uri=f"s3://{cxg_bucket.name}/{test_dataset_id}.cxg/",
            collection_visibility=CollectionVisibility.PUBLIC.name,
        )
        self.session.add(dataset)
        dataset = DbDataset(
            id="test_dataset_id_not_owner",
            revision=0,
            name="test_dataset_name_not_owner",
            organism={"ontology_term_id": "test_obo", "label": "test_organism"},
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=["test_sex", "test_sex2"],
            ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            collection_id="test_collection_id_not_owner",
            collection_visibility=CollectionVisibility.PRIVATE.name,
            explorer_url="test_url",
            explorer_s3_uri=f"s3://{cxg_bucket.name}/{test_dataset_id}.cxg/",
        )
        self.session.add(dataset)
        self.session.commit()

    def _create_test_dataset_artifacts(self):
        dataset_artifact = DbDatasetArtifact(
            id="test_dataset_artifact_id",
            dataset_id="test_dataset_id",
            filename="test_filename",
            filetype=DatasetArtifactFileType.H5AD.name,
            type=DatasetArtifactType.ORIGINAL.name,
            user_submitted=True,
            s3_uri=self.real_s3_file if self.real_data else self.fake_s3_file,
        )
        self.session.add(dataset_artifact)
        self.session.commit()

    def _create_test_dataset_processing_status(self):
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_id",
            dataset_id="test_dataset_id",
            processing_status=ProcessingStatus.PENDING,
            upload_status=UploadStatus.UPLOADING,
            upload_progress=4 / 9,
            validation_status=ValidationStatus.NA,
            conversion_loom_status=ConversionStatus.NA,
            conversion_rds_status=ConversionStatus.NA,
            conversion_cxg_status=ConversionStatus.NA,
            conversion_anndata_status=ConversionStatus.NA,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
