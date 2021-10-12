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
    "test": "http://localhost:5000",
}

AUDIENCE = {
    "prod": "cellxgene.cziscience.com/",
    "staging": "cellxgene.staging.single-cell.czi.technology/",
    "test": "cellxgene.dev.single-cell.czi.technology/",
    "dev": "cellxgene.dev.single-cell.czi.technology/",
}


class BaseFunctionalTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]

    def setUp(self):
        super().setUp()
        self.get_auth_token()
        self.api = API_URL.get(self.deployment_stage)
        self.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_collection_id = "DNE"
        self.bad_file_id = "DNE"

    @classmethod
    def get_auth_token(cls):
        config = CorporaAuthConfig()

        response = requests.post(
            "https://czi-cellxgene-dev.us.auth0.com/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="password",
                username=config.test_account_username,
                password=config.test_account_password,
                audience=AUDIENCE.get(cls.deployment_stage),
                scope="openid profile email offline",
                client_id=config.client_id,
                client_secret=config.client_secret,
            ),
        )

        access_token = response.json()["access_token"]
        id_token = response.json()["id_token"]
        token = {"access_token": access_token, "id_token": id_token}
        cls.cookie = base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()

    def upload_and_wait(self, collection_uuid, dropbox_url, existing_dataset_id=None):
        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
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
                if cxg_status == rds_status == h5ad_status == "CONVERTED":
                    keep_trying = False
            if time.time() >= timer + 300:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 5 min. Check logs for dataset: {dataset_uuid}"
                )
            time.sleep(10)
        return dataset_uuid
