import base64
from dataclasses import dataclass
import json
import os
import time
import typing
import unittest

from unittest.mock import Mock, patch


from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.layers.api.portal_api import PortalApi
from backend.layers.api.router import portal_api
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProviderInterface
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface
from backend.layers.common.entities import CollectionMetadata, CollectionVersion, CollectionVersionId, DatasetMetadata, DatasetStatusGeneric, Link, OntologyTermId
from tests.unit.backend.corpora.api_server.mock_auth import MockOauthServer
from tests.unit.backend.corpora.api_server.config import TOKEN_EXPIRES
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock

from typing import List, Optional

import copy

@dataclass
class DatasetStatusUpdate:
    status_key: str
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
            "backend.corpora.lambdas.api.v1.authentication.assert_authorized_token",
            side_effect=mock_assert_authorized_token,
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


class NewBaseTest(BaseAuthAPITest):

    business_logic: BusinessLogic
    crossref_provider: CrossrefProviderInterface # Can be mocked from the tests
    uri_provider: UriProviderInterface

    sample_dataset_metadata: DatasetMetadata
    
    def setUp(self):
        super().setUp()
        os.environ.setdefault("APP_NAME", "corpora-api")

        database_provider = DatabaseProviderMock()
        self.crossref_provider = CrossrefProviderInterface()
        step_function_provider = StepFunctionProviderInterface()
        s3_provider = S3Provider()
        uri_provider = UriProviderInterface()
        uri_provider.validate = Mock(return_value=True) # By default, every link should be valid

        self.sample_dataset_metadata = DatasetMetadata(
            name = "test_dataset_name",
            organism = [OntologyTermId(label="test_organism_label", ontology_term_id="test_organism_term_id")],
            tissue = [OntologyTermId(label="test_tissue_label", ontology_term_id="test_tissue_term_id")],
            assay = [OntologyTermId(label="test_assay_label", ontology_term_id="test_assay_term_id")],
            disease = [OntologyTermId(label="test_disease_label", ontology_term_id="test_disease_term_id")],
            sex = [OntologyTermId(label="test_sex_label", ontology_term_id="test_sex_term_id")],
            self_reported_ethnicity = [OntologyTermId(label="test_self_reported_ethnicity_label", ontology_term_id="test_self_reported_ethnicity_term_id")],
            development_stage = [OntologyTermId(label="test_development_stage_label", ontology_term_id="test_development_stage_term_id")],
            cell_type = [OntologyTermId(label="test_cell_type_label", ontology_term_id="test_cell_type_term_id")],
            cell_count = 10,
            schema_version = "3.0.0",
            mean_genes_per_cell = 0.5,
            batch_condition = ["test_batch_1", "test_batch_2"],
            suspension_type = ["test_suspension_type"],
            donor_id = ["test_donor_1"],
            is_primary_data = "BOTH",
            x_approximate_distribution="normal",
        )

        self.business_logic = BusinessLogic(
            database_provider, 
            self.crossref_provider, 
            step_function_provider, 
            s3_provider, 
            uri_provider
        )

        pa = PortalApi(self.business_logic)

        import backend.layers.api.router
        backend.layers.api.router.portal_api = Mock(return_value=pa)

        from backend.corpora.api_server.app import app
        # from backend.corpora.api_server.app import configure_flask_app, create_flask_app
        # with EnvironmentSetup(dict(APP_NAME="corpora-api")):
        # flask_app = configure_flask_app(create_flask_app())
        # self.app = flask_app.test_client(use_cookies=False)
        self.app = app.test_client(use_cookies=False)

    def generate_unpublished_collection(self, owner = "test_user_id", links: List[Link] = [], add_datasets: int = 0) -> CollectionVersion:

        metadata = CollectionMetadata(
            "test_collection",
            "described",
            "john doe",
            "john.due@email.com", # typo is on purpose
            links
        )

        collection = self.business_logic.create_collection(owner, metadata)

        for _ in range(add_datasets):
            
            metadata = copy.deepcopy(self.sample_dataset_metadata)
            # TODO: generate a real dataset, with artifact and processing status
            dataset_version_id, _ = self.business_logic.ingest_dataset(collection.version_id, "http://fake.url", None)
            self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
            # TODO: set a proper dataset status

        return collection

    # Public collections need to have at least one dataset!
    def generate_published_collection(self, owner = "test_user_id", links: List[Link] = [], add_datasets: int = 1):
        unpublished_collection = self.generate_unpublished_collection(owner, links, add_datasets=add_datasets)
        self.business_logic.publish_collection_version(unpublished_collection.version_id)
        return self.business_logic.get_collection_version(unpublished_collection.version_id)

    def generate_dataset(
        self, 
        owner: str = "test_user_id", 
        # collection_version_id: Optional[CollectionVersionId] = None, # TODO: probably remove
        metadata: Optional[DatasetMetadata] = None, 
        statuses: List[DatasetStatusUpdate] = [], 
        artifacts: List[DatasetArtifactUpdate] = [],
        publish: bool = False,
    ):
        """
        Convenience method for generating a dataset. Also generates an unpublished collection if needed.
        """
        # if not collection_version_id:
        collection = self.generate_unpublished_collection(owner)
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, "http://fake.url", None)
        if not metadata:
            metadata = copy.deepcopy(self.sample_dataset_metadata)
        self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
        for status in statuses:
            self.business_logic.update_dataset_version_status(dataset_version_id, status.status_key, status.status)
        artifact_ids = []
        for artifact in artifacts:
            artifact_ids.append(self.business_logic.add_dataset_artifact(dataset_version_id, artifact.type, artifact.uri))
        if publish:
            self.business_logic.publish_collection_version(collection.version_id)
        explorer_url = f"http://base.url/{dataset_id}"
        return DatasetData(
            dataset_version_id.id, 
            dataset_id.id, 
            explorer_url, 
            collection.version_id.id, 
            collection.collection_id.id,
            [a.id for a in artifact_ids],
        )


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



class BaseAPITest(unittest.TestCase):
    """
    Provide access to the test APIs. All tests for APIs should inherit this class.
    """

    maxDiff = None  # Easier to compare json responses.

    def setUp(self):
        super().setUp()
        os.environ.setdefault("APP_NAME", "corpora")
        from backend.corpora.api_server.app import app
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @staticmethod
    def remove_timestamps(body: dict, remove: typing.List[str] = None) -> dict:
        """
        A helper function to remove timestamps from the response body.
        :param body: The decoded json response body
        :param remove: Additional attributes to remove.
        :return: The decode json response body with timestamps removed.
        """
        defaults = ["created_at", "updated_at"]
        remove_attributes = remove + defaults if remove else defaults

        def _remove_timestamps(jrb, removing):
            if not isinstance(jrb, dict):
                return
            for rm in removing:
                jrb.pop(rm, None)
            for value in jrb.values():
                if isinstance(value, dict):
                    _remove_timestamps(value, removing)
                elif isinstance(value, list):
                    for list_value in value:
                        _remove_timestamps(list_value, removing)
            return jrb

        return _remove_timestamps(body, remove_attributes)


def mock_assert_authorized_token(token: str, audience: str = None):
    if token == "owner":
        return {"sub": "test_user_id", "email": "fake_user@email.com", "scope": []}
    elif token == "not_owner":
        return {"sub": "someone_else", "email": "fake_user@email.com", "scope": []}
    elif token == "super":
        return {"sub": "super", "email": "fake_user@email.com", "scope": ["write:collections"]}
    else:
        raise Exception()


def get_cxguser_token(user="owner"):
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



class AuthServerAPITest(BaseAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        (mock_oauth_server, auth_config) = cls.get_mock_server_and_auth_config()
        cls.mock_oauth_server = mock_oauth_server
        cls.auth_config = auth_config

    @staticmethod
    def get_mock_server_and_auth_config(additional_scope=None, token_duration=0):
        mock_oauth_server = MockOauthServer(additional_scope, token_duration)
        mock_oauth_server.start()
        assert mock_oauth_server.server_okay

        # Use the CorporaAuthConfig used by the app
        auth_config = CorporaAuthConfig()

        os.environ["API_BASE_URL"] = f"http://localhost:{mock_oauth_server.port}"
        # Overwrite the environment's auth config with our oidc server's config.
        authconfig = {
            "api_base_url": f"http://localhost:{mock_oauth_server.port}",
            "callback_base_url": auth_config.callback_base_url,
            "redirect_to_frontend": auth_config.redirect_to_frontend,
            "client_id": auth_config.client_id,
            "client_secret": auth_config.client_secret,
            "audience": auth_config.audience,
            "api_audience": auth_config.api_audience,
            "cookie_name": auth_config.cookie_name,
            "auth0_domain": f"localhost:{mock_oauth_server.port}",
            "curation_audience": auth_config.audience,
        }
        auth_config.set(authconfig)
        return (mock_oauth_server, auth_config)

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()
