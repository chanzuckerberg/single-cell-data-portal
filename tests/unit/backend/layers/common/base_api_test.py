import base64
import copy
import json
import os
import time
import typing
import unittest
from dataclasses import dataclass
from typing import List, Optional
from unittest.mock import Mock, patch

from backend.layers.api.portal_api import PortalApi
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionMetadata,
    CollectionVersion,
    DatasetMetadata,
    DatasetStatusGeneric,
    DatasetStatusKey,
    Link,
    OntologyTermId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.persistence.persistence_mock import DatabaseProviderMock
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProviderInterface
from tests.unit.backend.api_server.config import TOKEN_EXPIRES


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


class BaseAuthAPITest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.mock_assert_authorized_token = patch(
            "backend.portal.api.app.v1.authentication.assert_authorized_token",
            side_effect=self._mock_assert_authorized_token,
        )
        self.mock_assert_authorized_token.start()

    def tearDown(self):
        super().tearDown()
        self.mock_assert_authorized_token.stop()

    def make_owner_header(self):
        return {"Authorization": "Bearer " + "owner", "Content-Type": "application/json"}

    def make_super_curator_header(self):
        return {"Authorization": "Bearer " + "super", "Content-Type": "application/json"}

    def make_not_owner_header(self):
        return {"Authorization": "Bearer " + "not_owner", "Content-Type": "application/json"}

    def _mock_assert_authorized_token(self, token: str, audience: str = None):
        if token == "owner":
            return {"sub": "test_user_id", "email": "fake_user@email.com", "scope": []}
        elif token == "not_owner":
            return {"sub": "someone_else", "email": "fake_user@email.com", "scope": []}
        elif token == "super":
            return {"sub": "super", "email": "fake_user@email.com", "scope": ["write:collections"]}
        else:
            raise Exception()


class NewBaseTest(BaseAuthAPITest):

    business_logic: BusinessLogic
    crossref_provider: CrossrefProviderInterface  # Can be mocked from the tests
    uri_provider: UriProviderInterface

    sample_dataset_metadata: DatasetMetadata

    def setUp(self):
        super().setUp()
        os.environ.setdefault("APP_NAME", "corpora-api")

        database_provider = DatabaseProviderMock()
        # uncomment below and comment out above to run as integration tests, alongside a docker postgres instance
        # database_provider = DatabaseProvider()
        # database_provider._drop()
        # database_provider._create()

        self.crossref_provider = CrossrefProviderInterface()
        step_function_provider = StepFunctionProviderInterface()
        self.s3_provider = S3Provider()
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

        self.business_logic = BusinessLogic(
            database_provider, self.crossref_provider, step_function_provider, self.s3_provider, self.uri_provider
        )

        pa = PortalApi(self.business_logic)

        import backend.layers.api.router

        backend.layers.api.router.portal_api = Mock(return_value=pa)

        from backend.api_server.app import app

        # from backend.corpora.api_server.app import configure_flask_app, create_flask_app
        # with EnvironmentSetup(dict(APP_NAME="corpora-api")):
        # flask_app = configure_flask_app(create_flask_app())
        # self.app = flask_app.test_client(use_cookies=False)
        self.app = app.test_client(use_cookies=False)

    def generate_unpublished_collection(
        self, owner="test_user_id", links: List[Link] = [], add_datasets: int = 0
    ) -> CollectionVersion:

        metadata = CollectionMetadata("test_collection", "described", "john doe", "john.doe@email.com", links)

        collection = self.business_logic.create_collection(owner, metadata)

        for _ in range(add_datasets):

            metadata = copy.deepcopy(self.sample_dataset_metadata)
            # TODO: generate a real dataset, with artifact and processing status
            dataset_version_id, _ = self.business_logic.ingest_dataset(collection.version_id, "http://fake.url", None)
            self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
            # TODO: set a proper dataset status

        return self.business_logic.get_collection_version(collection.version_id)

    # Public collections need to have at least one dataset!
    def generate_published_collection(
            self, owner="test_user_id", links: List[Link] = [], add_datasets: int = 1
    ) -> CollectionVersion:
        unpublished_collection = self.generate_unpublished_collection(owner, links, add_datasets=add_datasets)
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
        statuses: List[DatasetStatusUpdate] = [], 
        artifacts: List[DatasetArtifactUpdate] = [],
        publish: bool = False,
    ) -> DatasetData:
        """
        Convenience method for generating a dataset. Also generates an unpublished collection if needed.
        """
        if not collection_version:
            collection_version = self.generate_unpublished_collection(owner)
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection_version.version_id, "http://fake.url", None)
        if not metadata:
            metadata = copy.deepcopy(self.sample_dataset_metadata)
        self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
        for status in statuses:
            self.business_logic.update_dataset_version_status(dataset_version_id, status.status_key, status.status)
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

    def get_cxguser_token(self, user="owner"):
        """
        Generated an auth token for testing.
        :param user: the type of use the token will simulate.
        :return:
        """
        cxguser = base64.b64encode(
            json.dumps(
                {
                    "access_token": user,
                    "refresh_token": f"random-{time.time()}",
                    "scope": "openid profile email offline",
                    "expires_in": TOKEN_EXPIRES,
                    "token_type": "Bearer",
                    "expires_at": TOKEN_EXPIRES,
                }
            ).encode("utf8")
        ).decode("utf8")
        return f"cxguser={cxguser}"
