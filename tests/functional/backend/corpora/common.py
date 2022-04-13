import unittest
import os
import requests
import base64
import json
import time

from backend.corpora.common.corpora_config import CorporaAuthConfig

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
    "test": "http://localhost:5005",
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
        token = cls.get_auth_token(cls.config.test_account_username, cls.config.test_account_password)
        cls.curator_cookie = cls.make_cookie(token)

    def setUp(self):
        super().setUp()
        self.api = API_URL.get(self.deployment_stage)
        self.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_collection_id = "DNE"
        self.bad_file_id = "DNE"

    @classmethod
    def get_auth_token(cls, username: str, password: str, additional_claims: list = None):
        standard_claims = "openid profile email offline"
        if additional_claims:
            additional_claims.append(standard_claims)
            claims = " ".join(additional_claims)
        else:
            claims = standard_claims
        response = requests.post(
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
        access_token = response.json()["access_token"]
        id_token = response.json()["id_token"]
        token = {"access_token": access_token, "id_token": id_token}
        return token

    @staticmethod
    def make_cookie(token: dict) -> str:
        return base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()

    def upload_and_wait(self, collection_uuid, dropbox_url, existing_dataset_id=None):
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        body = {"url": dropbox_url}

        if existing_dataset_id is None:
            res = requests.post(
                f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
            )
        else:
            body["id"] = existing_dataset_id
            res = requests.put(
                f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
            )

        res.raise_for_status()
        dataset_uuid = json.loads(res.content)["dataset_uuid"]
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/datasets/{dataset_uuid}", headers=headers)

        keep_trying = True
        timer = time.time()
        while keep_trying:
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
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
                    f"Dataset upload or conversion timed out after 10 min. Check logs for dataset: {dataset_uuid}"
                )
            time.sleep(10)
        return dataset_uuid
