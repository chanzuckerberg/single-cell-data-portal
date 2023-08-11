import base64
import json
import os
import time
import unittest

import requests  # type: ignore
from requests.adapters import HTTPAdapter, Response  # type: ignore
from requests.packages.urllib3.util import Retry  # type: ignore

from backend.common.corpora_config import CorporaAuthConfig

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
    "test": "https://localhost:5000",
}

AUDIENCE = {
    "prod": "api.cellxgene.cziscience.com",
    "staging": "api.cellxgene.staging.single-cell.czi.technology",
    "test": "api.cellxgene.dev.single-cell.czi.technology",
    "dev": "api.cellxgene.dev.single-cell.czi.technology",
}


class BaseFunctionalTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]
        cls.config = CorporaAuthConfig()
        cls.session = requests.Session()
        # apply retry config to idempotent http methods we use + POST requests, which are currently all either
        # idempotent (wmg queries) or low risk to rerun in dev/staging. Update if this changes in functional tests.
        retry_config = Retry(
            total=7,
            backoff_factor=2,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods={"DELETE", "GET", "HEAD", "PUT" "POST"},
        )
        cls.session.mount("https://", HTTPAdapter(max_retries=retry_config))
        token = cls.get_auth_token(cls.config.test_account_username, cls.config.test_account_password)
        cls.curator_cookie = cls.make_cookie(token)
        cls.api = API_URL.get(cls.deployment_stage)
        cls.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        cls.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        cls.bad_collection_id = "DNE"
        cls.bad_file_id = "DNE"
        cls.curation_api_access_token = cls.get_curation_api_access_token()

    @classmethod
    def get_curation_api_access_token(cls):
        response = cls.session.post(
            f"https://api.cellxgene.{cls.deployment_stage}.single-cell.czi.technology/curation/v1/auth/token",
            headers={"x-api-key": cls.config.super_curator_api_key},
        )
        return response.json()["access_token"]

    @classmethod
    def get_auth_token(cls, username: str, password: str, additional_claims: list = None):  # type: ignore
        standard_claims = "openid profile email offline"
        if additional_claims:
            additional_claims.append(standard_claims)
            claims = " ".join(additional_claims)
        else:
            claims = standard_claims
        response = cls.session.post(  # type: ignore
            "https://czi-cellxgene-dev.us.auth0.com/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="password",
                username=username,
                password=password,
                audience=AUDIENCE.get(cls.deployment_stage),  # type: ignore
                scope=claims,
                client_id=cls.config.client_id,  # type: ignore
                client_secret=cls.config.client_secret,  # type: ignore
            ),
        )
        access_token = response.json()["access_token"]
        id_token = response.json()["id_token"]
        token = {"access_token": access_token, "id_token": id_token}
        return token

    @staticmethod
    def make_cookie(token: dict) -> str:
        return base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()

    def upload_and_wait(self, collection_id, dropbox_url, existing_dataset_id=None, cleanup=True):
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        body = {"url": dropbox_url}

        if existing_dataset_id is None:
            res = self.session.post(
                f"{self.api}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=headers
            )
        else:
            body["id"] = existing_dataset_id
            res = self.session.put(
                f"{self.api}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=headers
            )

        res.raise_for_status()
        dataset_id = json.loads(res.content)["dataset_id"]

        if cleanup:
            self.addCleanup(requests.delete, f"{self.api}/dp/v1/datasets/{dataset_id}", headers=headers)

        keep_trying = True
        timer = time.time()
        while keep_trying:
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=headers)
            res.raise_for_status()
            data = json.loads(res.content)
            upload_status = data["upload_status"]
            if upload_status == "UPLOADED":
                cxg_status = data.get("cxg_status")
                rds_status = data.get("rds_status")
                h5ad_status = data.get("h5ad_status")
                processing_status = data.get("processing_status")
                if cxg_status == rds_status == h5ad_status == "UPLOADED" and processing_status == "SUCCESS":
                    keep_trying = False
            if time.time() >= timer + 600:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 10 min. Check logs for dataset: {dataset_id}"
                )
            time.sleep(10)
        return dataset_id

    def assertStatusCode(self, actual: int, expected_response: Response):
        request_id = expected_response.headers.get("X-Request-Id")
        self.assertEqual(actual, expected_response.status_code, msg=f"{request_id=}")
