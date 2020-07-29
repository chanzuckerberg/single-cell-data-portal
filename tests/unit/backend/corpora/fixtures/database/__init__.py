import os
import random
import string
import sys

from sqlalchemy import func

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.scripts.create_db import create_db
from backend.corpora.common.utils.db_utils import DbUtils
from backend.corpora.common.corpora_orm import (
    ProjectStatus,
    ProcessingState,
    ValidationState,
    ProjectLinkType,
    DatasetArtifactType,
    DatasetArtifactFileType,
    DbProject,
    DbProjectLink,
    DbDataset,
    DbDatasetArtifact,
    DbContributor,
    DbDatasetContributor,
    DbUser,
    DbDeploymentDirectory,
    Base
)
from server.db.cellxgene_orm import CellxGeneUser, CellxGeneDataset, Annotation
from server.db.cellxgene_orm import Base as cxgBase


class TestDatabase:
    def __init__(self):
        create_db()
        self.db = DbUtils()
        self._populate_test_data()
        del self.db

    def _populate_test_data(self):
        self._create_test_users()
        self._create_test_projects()
        self._create_test_project_links()
        self._create_test_datasets()
        self._create_test_dataset_artifacts()
        self._create_test_contributors()
        self._create_test_cellxgene_user()
        self._create_test_cellxgene_dataset()
        self._create_test_annotation()
        self._populate_test_data_many()

    def _create_test_users(self):
        user = DbUser(id="test_user_id", name="test_user", email="test_email")
        self.db.session.add(user)
        self.db.session.commit()

    def _create_test_projects(self):
        project = DbProject(
            id="test_project_id",
            status=ProjectStatus.LIVE.name,
            owner="test_user_id",
            name="test_project",
            description="test_description",
            s3_bucket="test_s3_bucket",
            tc_uri="test_tc_uri",
            needs_attestation=False,
            processing_state=ProcessingState.NA.name,
            validation_state=ValidationState.NOT_VALIDATED.name,
        )
        self.db.session.add(project)
        self.db.session.commit()

    def _create_test_project_links(self):
        project_link = DbProjectLink(
            id="test_project_link_id",
            project_id="test_project_id",
            project_status=ProjectStatus.LIVE.name,
            link_name="test_link_name",
            link_url="test_url",
            link_type=ProjectLinkType.RAW_DATA.name,
        )
        self.db.session.add(project_link)
        self.db.session.commit()

    def _create_test_datasets(self):
        test_dataset_id = "test_dataset_id"
        dataset = DbDataset(
            id=test_dataset_id,
            revision=0,
            name="test_dataset_name",
            organism="test_organism",
            organism_ontology="test_organism_ontology",
            tissue="test_tissue",
            tissue_ontology="test_tissue_ontology",
            assay="test_assay",
            assay_ontology="test_assay_ontology",
            disease="test_disease",
            disease_ontology="test_disease_ontology",
            sex="test_sex",
            ethnicity="test_ethnicity",
            ethnicity_ontology="test_ethnicity_ontology",
            development_stage="test_development_stage",
            development_stage_ontology="test_development_stage_ontology",
            source_data_location="test_source_data_location",
            preprint_doi="test_preprint_doi",
            publication_doi="test_publication_doi",
            project_id="test_project_id",
            project_status=ProjectStatus.LIVE.name,
        )
        self.db.session.add(dataset)
        self.db.session.commit()

        deployment_directory = DbDeploymentDirectory(
            id="test_deployment_directory_id", dataset_id=test_dataset_id, environment="test", url="test_url"
        )
        self.db.session.add(deployment_directory)
        self.db.session.commit()

    def _create_test_dataset_artifacts(self):
        dataset_artifact = DbDatasetArtifact(
            id="test_dataset_artifact_id",
            dataset_id="test_dataset_id",
            filename="test_filename",
            filetype=DatasetArtifactFileType.H5AD.name,
            type=DatasetArtifactType.ORIGINAL.name,
            user_submitted=True,
            s3_uri="test_s3_uri",
        )
        self.db.session.add(dataset_artifact)
        self.db.session.commit()

    def _create_test_contributors(self):
        contributor = DbContributor(
            id="test_contributor_id", name="test_contributor_name", institution="test_institution", email="test_email"
        )
        self.db.session.add(contributor)
        self.db.session.commit()

        dataset_contributor = DbDatasetContributor(
            id="test_dataset_contributor_id", contributor_id="test_contributor_id", dataset_id="test_dataset_id"
        )
        self.db.session.add(dataset_contributor)
        self.db.session.commit()

    def _populate_test_data_many(self):
        self._create_test_cellxgene_users()
        self._create_test_cellxgene_datasets()
        self._create_test_annotations()

    def _create_test_cellxgene_user(self):
        user = CellxGeneUser(id="test_user_id")
        self.db.session.add(user)
        self.db.session.commit()

    def _create_test_cellxgene_dataset(self):
        dataset = CellxGeneDataset(
            id="test_dataset_id",
            name="test_dataset",
        )
        self.db.session.add(dataset)
        self.db.session.commit()

    def _create_test_annotation(self):
        annotation = Annotation(
            id="test_annotation_id",
            tiledb_uri="tiledb_uri",
            user_id="test_user_id",
            dataset_id="test_dataset_id"
        )
        self.db.session.add(annotation)
        self.db.session.commit()

    @staticmethod
    def get_random_string():
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(12))

    def _create_test_cellxgene_users(self, user_count: int = 10):
        users = []
        for i in range(user_count):
            users.append(CellxGeneUser(id=self.get_random_string()))
        self.db.session.add_all(users)
        self.db.session.commit()

    def _create_test_cellxgene_datasets(self, dataset_count: int = 10):
        datasets = []
        for i in range(dataset_count):
            datasets.append(CellxGeneDataset(id=self.get_random_string(), name=self.get_random_string()))
        self.db.session.add_all(datasets)
        self.db.session.commit()

    def _create_test_annotations(self, annotation_count: int = 10):
        annotations = []
        for i in range(annotation_count):
            dataset = self.order_by_random(CellxGeneDataset)
            user = self.order_by_random(DbUser)
            annotations.append(Annotation(
                id=self.get_random_string(),
                tiledb_uri=self.get_random_string(),
                user_id=user.id,
                dataset_id=dataset.id
            ))
        self.db.session.add_all(annotations)
        self.db.session.commit()

    def order_by_random(self, table: cxgBase):
        return self.db.session.query(table).order_by(func.random()).first()
