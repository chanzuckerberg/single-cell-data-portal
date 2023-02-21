import numpy as np

from backend.common.corpora_orm import (
    CollectionLinkType,
    CollectionVisibility,
    ConversionStatus,
    DatasetArtifactFileType,
    DbCollection,
    DbCollectionLink,
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    DbGeneset,
    IsPrimaryData,
    ProcessingStatus,
    UploadStatus,
    ValidationStatus,
    XApproximateDistribution,
)
from backend.common.utils.db_session import DBSessionMaker
from backend.scripts.create_db import create_db
from tests.unit.backend.fixtures import config


class DatabaseManager:
    is_initialized = False

    @classmethod
    def initialize_db(cls):
        if cls.is_initialized:
            return
        testdb = DatabaseFixture()
        testdb.create_db()
        testdb.populate_test_data()
        cls.is_initialized = True


class DatabaseFixture:
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
        # TODO update for the redesign
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
            id="test_collection_id_revision",
            visibility=CollectionVisibility.PRIVATE.name,
            revision_of="test_collection_id",
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
            id="test_collection_id_public_for_revision_one",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="User1",
            name="test_collection_id_public_for_revision_one",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
        )
        self.session.add(collection)
        collection = DbCollection(
            id="test_collection_id_public_for_revision_two",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="User1",
            name="test_collection_id_public_for_revision_two",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
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
        collection = DbCollection(
            id="test_collection_with_link",
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
            id="test_collection_with_link_and_dataset_changes",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
        )
        self.session.add(collection)
        self.session.commit()

    def _create_test_geneset(self):
        geneset = DbGeneset(
            id="test_geneset",
            name="test_geneset",
            description="this is a geneset",
            collection_id="test_collection_id",
        )

        self.session.add(geneset)
        dataset = self.session.query(DbDataset).get("test_dataset_id")
        geneset = DbGeneset(
            id="test_geneset_with_dataset",
            name="test_geneset_with_dataset",
            description="this is a geneset with a dataset",
            collection_id="test_collection_id",
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
                    link_name=f"test_{link_type.value}_link_name",
                    link_url=f"http://test_{link_type.value}_url.place",
                    link_type=link_type.name,
                )
            )
            self.session.add(
                DbCollectionLink(
                    id=f"test_collection_no_name_{link_type.value}_link_id",
                    collection_id="test_collection_id",
                    link_url=f"http://test_no_link_name_{link_type.value}_url.place",
                    link_type=link_type.name,
                )
            )
            self.session.add(
                DbCollectionLink(
                    id=f"test_publish_revision_with_link__{link_type.value}_link",
                    collection_id="test_collection_with_link",
                    link_name=f"test_{link_type.value}_link_name",
                    link_url=f"http://test_link_{link_type.value}_url.place",
                    link_type=link_type.name,
                )
            )
            self.session.add(
                DbCollectionLink(
                    id=f"test_publish_revision_with_link_and_dataset_changes__{link_type.value}_link",
                    collection_id="test_collection_with_link_and_dataset_changes",
                    link_name=f"test_{link_type.value}_link_name",
                    link_url=f"http://test_link_{link_type.value}_url.place",
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
            organism=[{"ontology_term_id": "test_obo", "label": "test_organism"}],
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=[
                {"label": "test_sex", "ontology_term_id": "test_obo"},
                {"label": "test_sex2", "ontology_term_id": "test_obp"},
            ],
            self_reported_ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            collection_id="test_collection_id",
            explorer_url="test_url",
            cell_type=[{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
            is_primary_data=IsPrimaryData.PRIMARY.name,
            x_approximate_distribution=XApproximateDistribution.NORMAL.name,
            batch_condition=np.array(["batchA", "batchB"], dtype="object"),
            donor_id=["donor_1", "donor_2"],
            suspension_type=["nucleus"],
            schema_version="3.0.0",
        )
        self.session.add(dataset)
        dataset = DbDataset(
            id="test_dataset_for_revisions_one",
            revision=0,
            name="test_dataset_name",
            organism=[{"ontology_term_id": "test_obo", "label": "test_organism"}],
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=[
                {"label": "test_sex", "ontology_term_id": "test_obo"},
                {"label": "test_sex2", "ontology_term_id": "test_obp"},
            ],
            self_reported_ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            donor_id=["donor_1", "donor_2"],
            suspension_type=["nucleus"],
            collection_id="test_collection_id_public_for_revision_one",
            explorer_url="test_url",
            schema_version="3.0.0",
        )
        self.session.add(dataset)
        dataset = DbDataset(
            id="test_dataset_for_revisions_two",
            revision=0,
            name="test_dataset_name",
            organism=[{"ontology_term_id": "test_obo", "label": "test_organism"}],
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=[
                {"label": "test_sex", "ontology_term_id": "test_obo"},
                {"label": "test_sex2", "ontology_term_id": "test_obp"},
            ],
            self_reported_ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            donor_id=["donor_1", "donor_2"],
            suspension_type=["nucleus"],
            collection_id="test_collection_id_public_for_revision_two",
            explorer_url="test_url",
            schema_version="3.0.0",
        )
        self.session.add(dataset)
        dataset = DbDataset(
            id="test_dataset_id_not_owner",
            revision=0,
            name="test_dataset_name_not_owner",
            organism=[{"ontology_term_id": "test_obo", "label": "test_organism"}],
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=[
                {"label": "test_sex", "ontology_term_id": "test_obo"},
                {"label": "test_sex2", "ontology_term_id": "test_obp"},
            ],
            self_reported_ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            donor_id=["donor_1", "donor_2"],
            suspension_type=["nucleus"],
            collection_id="test_collection_id_not_owner",
            explorer_url="test_url",
            cell_type=[{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
            is_primary_data=IsPrimaryData.PRIMARY.name,
            x_approximate_distribution=XApproximateDistribution.NORMAL.name,
            schema_version="3.0.0",
        )
        self.session.add(dataset)
        dataset = DbDataset(
            id="test_publish_revision_with_links__revision_dataset",
            original_id=test_dataset_id,
            revision=0,
            name="test_dataset_name_revised",
            organism=[{"ontology_term_id": "test_obo", "label": "test_organism"}],
            tissue=[{"ontology_term_id": "test_obo", "label": "test_tissue"}],
            assay=[{"ontology_term_id": "test_obo", "label": "test_assay"}],
            disease=[
                {"ontology_term_id": "test_obo", "label": "test_disease"},
                {"ontology_term_id": "test_obp", "label": "test_disease2"},
                {"ontology_term_id": "test_obq", "label": "test_disease3"},
            ],
            sex=[
                {"label": "test_sex", "ontology_term_id": "test_obo"},
                {"label": "test_sex2", "ontology_term_id": "test_obp"},
            ],
            self_reported_ethnicity=[{"ontology_term_id": "test_obo", "label": "test_ethnicity"}],
            development_stage=[{"ontology_term_id": "test_obo", "label": "test_development_stage"}],
            donor_id=["donor_1", "donor_2"],
            suspension_type=["nucleus"],
            collection_id="test_collection_id_revision",
            explorer_url="test_url_revised",
            cell_type=[{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
            is_primary_data=IsPrimaryData.PRIMARY.name,
            x_approximate_distribution=XApproximateDistribution.NORMAL.name,
            schema_version="3.0.0",
        )
        self.session.add(dataset)
        self.session.commit()

    def _create_test_dataset_artifacts(self):
        dataset_artifact = DbDatasetArtifact(
            id="test_dataset_artifact_id",
            dataset_id="test_dataset_id",
            filename="test_filename",
            filetype=DatasetArtifactFileType.H5AD.name,
            user_submitted=True,
            s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
        )
        self.session.add(dataset_artifact)
        self.session.commit()

        dataset_artifact = DbDatasetArtifact(
            id="test_dataset_artifact_id_cxg",
            dataset_id="test_dataset_id",
            filename="test_filename",
            filetype=DatasetArtifactFileType.CXG.name,
            user_submitted=True,
            s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
        )
        self.session.add(dataset_artifact)
        self.session.commit()

        # Revision 1

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_one_CXG_id",
                dataset_id="test_dataset_for_revisions_one",
                filename="test_filename",
                filetype=DatasetArtifactFileType.CXG.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_one_H5AD_id",
                dataset_id="test_dataset_for_revisions_one",
                filename="test_filename",
                filetype=DatasetArtifactFileType.H5AD.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        self.session.commit()

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_one_RDS_id",
                dataset_id="test_dataset_for_revisions_one",
                filename="test_filename",
                filetype=DatasetArtifactFileType.RDS.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        # Revision 2

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_two_CXG_id",
                dataset_id="test_dataset_for_revisions_two",
                filename="test_filename",
                filetype=DatasetArtifactFileType.CXG.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_two_H5AD_id",
                dataset_id="test_dataset_for_revisions_two",
                filename="test_filename",
                filetype=DatasetArtifactFileType.H5AD.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        self.session.add(
            DbDatasetArtifact(
                id="test_dataset_artifact_for_revision_two_RDS_id",
                dataset_id="test_dataset_for_revisions_two",
                filename="test_filename",
                filetype=DatasetArtifactFileType.RDS.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

        self.session.add(
            DbDatasetArtifact(
                id="test_publish_revision_with_links__revision_artifact",
                dataset_id="test_publish_revision_with_links__revision_dataset",
                filename="test_filename",
                filetype=DatasetArtifactFileType.H5AD.name,
                user_submitted=True,
                s3_uri=config.real_s3_file if self.real_data else config.fake_s3_file,
            )
        )

        self.session.commit()

    def _create_test_dataset_processing_status(self):
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_id",
            dataset_id="test_dataset_id",
            processing_status=ProcessingStatus.PENDING,
            upload_status=UploadStatus.UPLOADING,
            upload_progress=4 / 9,
            validation_status=ValidationStatus.NA,
            rds_status=ConversionStatus.NA,
            cxg_status=ConversionStatus.NA,
            h5ad_status=ConversionStatus.NA,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_for_revisions_one_id",
            dataset_id="test_dataset_for_revisions_one",
            processing_status=ProcessingStatus.SUCCESS,
            upload_status=UploadStatus.UPLOADED,
            upload_progress=1,
            validation_status=ValidationStatus.VALID,
            rds_status=ConversionStatus.CONVERTED,
            cxg_status=ConversionStatus.CONVERTED,
            h5ad_status=ConversionStatus.CONVERTED,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_for_revisions_two_id",
            dataset_id="test_dataset_for_revisions_two",
            processing_status=ProcessingStatus.SUCCESS,
            upload_status=UploadStatus.UPLOADED,
            upload_progress=1,
            validation_status=ValidationStatus.VALID,
            rds_status=ConversionStatus.CONVERTED,
            cxg_status=ConversionStatus.CONVERTED,
            h5ad_status=ConversionStatus.CONVERTED,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_for_revision_with_links",
            dataset_id="test_publish_revision_with_links__revision_dataset",
            processing_status=ProcessingStatus.SUCCESS,
            upload_status=UploadStatus.UPLOADED,
            upload_progress=1,
            validation_status=ValidationStatus.VALID,
            rds_status=ConversionStatus.CONVERTED,
            cxg_status=ConversionStatus.CONVERTED,
            h5ad_status=ConversionStatus.CONVERTED,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
        dataset_processing_status = DbDatasetProcessingStatus(
            id="test_dataset_processing_status_for_not_owner",
            dataset_id="test_dataset_id_not_owner",
            processing_status=ProcessingStatus.PENDING,
            upload_status=UploadStatus.UPLOADED,
            upload_progress=1,
            validation_status=ValidationStatus.VALID,
            rds_status=ConversionStatus.CONVERTED,
            cxg_status=ConversionStatus.CONVERTED,
            h5ad_status=ConversionStatus.CONVERTED,
        )
        self.session.add(dataset_processing_status)
        self.session.commit()
