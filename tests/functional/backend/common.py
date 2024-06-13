import base64
import json
import os
import time
import unittest
from typing import Optional

import requests
from requests.adapters import HTTPAdapter, Response, Retry

from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL, AUDIENCE


class BaseFunctionalTestCase(unittest.TestCase):
    session: requests.Session
    config: CorporaAuthConfig
    deployment_stage: str

    @classmethod
    def setUpClass(cls, smoke_tests: bool = False):
        super().setUpClass()
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]
        cls.config = CorporaAuthConfig()
        cls.test_dataset_uri = (
            "https://www.dropbox.com/scl/fi/y50umqlcrbz21a6jgu99z/"
            "5_0_0_example_valid.h5ad?rlkey=s7p6ybyx082hswix26hbl11pm&dl=0"
        )
        cls.test_visium_dataset_uri = (
            "https://www.dropbox.com/scl/fi/lmhue0va6ihk50ivp26da/visium_small.h5ad?rlkey=n0f"
            "o4dyi1ah7ckg9kgzwlhm8s&st=p7jzej8j&dl=0"
        )
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
        if cls.deployment_stage == "rdev":
            cls.get_oauth2_proxy_access_token()
        if smoke_tests:
            username, password = cls.config.test_account_username, cls.config.test_account_password
        else:
            username, password = cls.config.functest_account_username, cls.config.functest_account_password
        token = cls.get_auth_token(username, password)
        cls.curator_cookie = cls.make_cookie(token)
        cls.api = API_URL.get(cls.deployment_stage)
        cls.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        cls.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        cls.bad_collection_id = "DNE"
        cls.bad_file_id = "DNE"
        cls.curation_api_access_token = cls.get_curation_api_access_token()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.session.close()

    @classmethod
    def get_curation_api_access_token(cls):
        response = cls.session.post(
            f"{cls.api}/curation/v1/auth/token",
            headers={"x-api-key": cls.config.super_curator_api_key},
        )
        response.raise_for_status()
        return response.json()["access_token"]

    @classmethod
    def get_oauth2_proxy_access_token(cls):
        payload = {
            "client_id": cls.config.test_app_id,
            "client_secret": cls.config.test_app_secret,
            "grant_type": "client_credentials",
            "audience": "https://api.cellxgene.dev.single-cell.czi.technology/dp/v1/curator",
        }
        headers = {"content-type": "application/json"}

        res = cls.session.post("https://czi-cellxgene-dev.us.auth0.com/oauth/token", json=payload, headers=headers)
        res.raise_for_status()
        cls.proxy_access_token = res.json()["access_token"]
        cls.session.headers["Authorization"] = f"Bearer {cls.proxy_access_token}"

    @classmethod
    def get_auth_token(cls, username: str, password: str, additional_claims: Optional[list] = None):
        standard_claims = "openid profile email offline"
        if additional_claims:
            additional_claims.append(standard_claims)
            claims = " ".join(additional_claims)
        else:
            claims = standard_claims
        response = cls.session.post(
            "https://czi-cellxgene-dev.us.auth0.com/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="password",
                username=username,
                password=password,
                audience=AUDIENCE.get(cls.deployment_stage),
                scope=claims,
                client_id=cls.config.client_id,
                client_secret=cls.config.client_secret,
            ),
        )
        response.raise_for_status()
        access_token = response.json()["access_token"]
        id_token = response.json()["id_token"]
        token = {"access_token": access_token, "id_token": id_token}
        return token

    @staticmethod
    def make_cookie(token: dict) -> str:
        return base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()

    def upload_and_wait(self, collection_id, dropbox_url, existing_dataset_id=None, cleanup=True):
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        if self.deployment_stage == "rdev":
            headers["Authorization"] = f"Bearer {self.proxy_access_token}"
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

        self.assertStatusCode(requests.codes.accepted, res)

        keep_trying = True
        expected_upload_statuses = ["WAITING", "UPLOADING", "UPLOADED"]
        # conversion statuses can be `None` when/if we hit the status endpoint too early after an upload
        expected_conversion_statuses = ["CONVERTING", "CONVERTED", "FAILED", "UPLOADING", "UPLOADED", "NA", None]
        timer = time.time()
        while keep_trying:
            res = self.session.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=headers)
            res.raise_for_status()
            data = json.loads(res.content)
            upload_status = data["upload_status"]
            if upload_status:
                self.assertIn(upload_status, expected_upload_statuses)

            # conversion statuses only returned once uploaded
            if upload_status == "UPLOADED":
                cxg_status = data.get("cxg_status")
                rds_status = data.get("rds_status")
                h5ad_status = data.get("h5ad_status")
                self.assertIn(data.get("cxg_status"), expected_conversion_statuses)
                if cxg_status == "FAILED":
                    self.fail(f"CXG CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if rds_status == "FAILED":
                    self.fail(f"RDS CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if h5ad_status == "FAILED":
                    self.fail(f"Anndata CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if cxg_status == rds_status == h5ad_status == "UPLOADED":
                    keep_trying = False
            if time.time() >= timer + 600:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 10 min. Check logs for dataset: {dataset_id}"
                )
            time.sleep(10)
        return dataset_id

    def assertStatusCode(self, actual: int, expected_response: Response):
        request_id = expected_response.headers.get("X-Request-Id")
        assert actual == expected_response.status_code, f"{request_id=}"
