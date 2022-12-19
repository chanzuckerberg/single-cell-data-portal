import copy
import os
import typing
import unittest
from dataclasses import dataclass
from typing import List, Optional
from unittest.mock import Mock

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionMetadata,
    CollectionVersion,
    DatasetMetadata,
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetValidationStatus,
    Link,
    OntologyTermId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.persistence.persistence_mock import DatabaseProviderMock
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProviderInterface


@dataclass
class DatasetStatusUpdate:
    status_key: DatasetStatusKey
    status: DatasetStatusGeneric


@dataclass
class DatasetArtifactUpdate:
    type: str
    uri: str


@dataclass
class DatasetData:
    """
    Convenience class that returns all the information required by the tests.
    The ids are already stringified for convenience, since the API layer
    will work with strings (at least for now)
    """

    dataset_version_id: str
    dataset_id: str
    explorer_url: str
    collection_version_id: str
    collection_id: str
    artifact_ids: List[str]


class BaseTest(unittest.TestCase):

    business_logic: BusinessLogic
    crossref_provider: CrossrefProviderInterface  # Can be mocked from the tests
    uri_provider: UriProviderInterface

    sample_dataset_metadata: DatasetMetadata
    sample_collection_metadata: CollectionMetadata

    @classmethod
    def setUpClass(cls) -> None:
        cls.run_as_integration = True if os.environ.get("INTEGRATION_TEST", "false").lower() == "true" else False
        if cls.run_as_integration:
            database_uri = os.environ.get("DB_URI", "postgresql://postgres:secret@localhost")
            cls.database_provider = DatabaseProvider(database_uri=database_uri)
            cls.database_provider._drop_schema("persistence_schema")

    def setUp(self):
        super().setUp()
        os.environ.setdefault("APP_NAME", "corpora-api")

        if self.run_as_integration:
            self.database_provider._create_schema("persistence_schema")
        else:
            self.database_provider = DatabaseProviderMock()

        self.crossref_provider = CrossrefProviderInterface()
        step_function_provider = StepFunctionProviderInterface()
        self.s3_provider = S3ProviderInterface()
        self.uri_provider = UriProviderInterface()
        self.uri_provider.validate = Mock(return_value=True)  # By default, every link should be valid
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

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

        self.sample_collection_metadata = CollectionMetadata(
            "test_collection", "described", "john doe", "john.doe@email.com", []
        )

        self.business_logic = BusinessLogic(
            self.database_provider, self.crossref_provider, step_function_provider, self.s3_provider, self.uri_provider
        )

    def tearDown(self):
        super().tearDown()
        if self.run_as_integration:
            self.database_provider._drop_schema("persistence_schema")

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.run_as_integration:
            cls.database_provider._engine.dispose()

    def generate_unpublished_collection(
        self,
        owner="test_user_id",
        curator_name="Test User",
        links: List[Link] = [],
        add_datasets: int = 0,
        metadata=None,
    ) -> CollectionVersion:
        links = links if links else []
        if not metadata:
            metadata = copy.deepcopy(self.sample_collection_metadata)
            metadata.links = links
        collection = self.business_logic.create_collection(owner, curator_name, metadata)

        for _ in range(add_datasets):

            metadata = copy.deepcopy(self.sample_dataset_metadata)
            # TODO: generate a real dataset, with artifact and processing status
            dataset_version_id, _ = self.business_logic.ingest_dataset(collection.version_id, "http://fake.url", None)
            self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
            # TODO: set a proper dataset status

        return self.business_logic.get_collection_version(collection.version_id)

    # Public collections need to have at least one dataset!
    # Public collections need to have at least one dataset!
    def generate_published_collection(
        self,
        owner="test_user_id",
        links: List[Link] = [],
        add_datasets: int = 1,
        curator_name: str = "Jane Smith",
        metadata=None,
    ) -> CollectionVersion:
        unpublished_collection = self.generate_unpublished_collection(
            owner, curator_name, links, add_datasets=add_datasets, metadata=metadata
        )
        self.business_logic.publish_collection_version(unpublished_collection.version_id)
        return self.business_logic.get_collection_version(unpublished_collection.version_id)

    def generate_revision(self, collection_id: str) -> CollectionVersion:
        revision = self.business_logic.create_collection_version(collection_id)
        return self.business_logic.get_collection_version(revision.version_id)

    def generate_dataset(
        self,
        owner: str = "test_user_id",
        collection_version: Optional[CollectionVersion] = None,
        metadata: Optional[DatasetMetadata] = None,
        name: Optional[str] = None,
        statuses: List[DatasetStatusUpdate] = [],
        validation_message: str = None,
        artifacts: List[DatasetArtifactUpdate] = [],
        publish: bool = False,
    ) -> DatasetData:
        """
        Convenience method for generating a dataset. Also generates an unpublished collection if needed.
        """
        if not collection_version:
            collection_version = self.generate_unpublished_collection(owner)
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection_version.version_id, "http://fake.url", None
        )
        if not metadata:
            metadata = copy.deepcopy(self.sample_dataset_metadata)
        if name is not None:
            metadata.name = name
        self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
        for status in statuses:
            self.business_logic.update_dataset_version_status(dataset_version_id, status.status_key, status.status)
        if validation_message:
            self.business_logic.update_dataset_version_status(
                dataset_version_id,
                DatasetStatusKey.VALIDATION,
                DatasetValidationStatus.INVALID,
                validation_message=validation_message,
            )
        artifact_ids = []
        for artifact in artifacts:
            artifact_ids.append(
                self.business_logic.add_dataset_artifact(dataset_version_id, artifact.type, artifact.uri)
            )
        if publish:
            self.business_logic.publish_collection_version(collection_version.version_id)
        explorer_url = f"http://base.url/{dataset_version_id.id}.cxg/"
        return DatasetData(
            dataset_version_id.id,
            dataset_id.id,
            explorer_url,
            collection_version.version_id.id,
            collection_version.collection_id.id,
            [a.id for a in artifact_ids],
        )

    def _generate_id(self):
        return DatabaseProviderMock._generate_id()

    def remove_timestamps(self, body: dict, remove: typing.List[str] = None) -> dict:
        # TODO: implement as needed
        return body

    def generate_collection(self, links: List[dict] = None, **kwargs) -> CollectionVersion:
        """Generated a collection
        Adding to for compatibility with old tests
        """
        links = links if links else []
        links = [api_link_dict_to_link_class(lk) for lk in links]
        visibility = kwargs.pop("visibility", "PUBLIC")
        if visibility == "PUBLIC":
            return self.generate_published_collection(links=links, **kwargs)
        else:
            return self.generate_unpublished_collection(links=links, **kwargs)

    def generate_collection_revision(self, owner="test_user_id") -> CollectionVersion:
        published_collection = self.generate_published_collection(owner)
        return self.business_logic.create_collection_version(published_collection.collection_id)


def link_class_to_api_link_dict(link: Link) -> dict:
    return {"link_name": link.name, "link_type": link.type, "link_url": link.uri}


def api_link_dict_to_link_class(link: dict) -> Link:
    return Link(link["link_name"], link["link_type"], link["link_url"])
