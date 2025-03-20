#!/usr/bin/env python
import json
import os
import threading

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL, ATAC_SEQ_MANIFEST, DATASET_MANIFEST, VISIUM_DATASET_MANIFEST
from tests.functional.backend.utils import (
    get_auth_token,
    get_curation_api_access_token,
    make_cookie,
    make_proxy_auth_token,
    make_session,
    upload_manifest_and_wait,
)

# Amount to reduce chance of collision where multiple test instances select the same collection to test against
NUM_TEST_DATASETS = 3
NUM_TEST_COLLECTIONS = 10
TEST_ACCT_CONTACT_NAME = "Smoke Test User"
VISIUM_ACCT_CONTACT_NAME = "Visium Test User"
ATAC_SEQ_ACCT_CONTACT_NAME = "ATAC Seq Test User"


class SmokeTestsInitializer:
    def __init__(self):
        self.deployment_stage = os.environ["DEPLOYMENT_STAGE"]
        self.config = CorporaAuthConfig()
        proxy_auth_token = make_proxy_auth_token(self.config, self.deployment_stage)
        self.session = make_session(proxy_auth_token)
        self.api = API_URL.get(self.deployment_stage)
        username, password = self.config.test_account_username, self.config.test_account_password
        auth_token = get_auth_token(username, password, self.session, self.config, self.deployment_stage)
        self.curator_cookie = make_cookie(auth_token)
        self.curation_api_access_token = get_curation_api_access_token(self.session, self.api, self.config)
        self.manifests = [DATASET_MANIFEST, VISIUM_DATASET_MANIFEST, ATAC_SEQ_MANIFEST]
        self.headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

    def get_collection_count(self, contact_name=TEST_ACCT_CONTACT_NAME, expected_count=NUM_TEST_COLLECTIONS):
        res = self.session.get(f"{self.api}/curation/v1/collections?visiblity=PUBLIC", headers=self.headers)
        res.raise_for_status()
        data = json.loads(res.content)
        num_collections = 0
        for collection in data:
            if collection["contact_name"] == contact_name:
                num_collections += 1
            if num_collections == expected_count:
                return num_collections
        return num_collections

    def start_upload_thread(self, collection_id, manifest) -> threading.Thread:
        _thread = threading.Thread(
            target=upload_manifest_and_wait,
            args=(
                self.session,
                self.api,
                self.curation_api_access_token,
                self.curator_cookie,
                collection_id,
                manifest,
            ),
        )
        _thread.start()
        return _thread

    def create_and_publish_collection(
        self,
        contact_name=TEST_ACCT_CONTACT_NAME,
        collection_name="test collection",
        manifest=DATASET_MANIFEST,
        num_datasets=NUM_TEST_DATASETS,
    ):
        collection_id = self.create_collection(contact_name, collection_name)
        _threads = []
        for _ in range(num_datasets):
            _threads.append(self.start_upload_thread(collection_id, manifest))
        for _thread in _threads:
            _thread.join()
        self.publish_collection(collection_id)
        print(f"created and published collection {collection_id}")

    def create_collection(self, contact_name, collection_name):
        data = {
            "contact_email": "example@gmail.com",
            "contact_name": contact_name,
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "http://protocol.com"}],
            "name": collection_name,
        }

        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=self.headers)
        res.raise_for_status()
        data = json.loads(res.content)
        return data["collection_id"]

    def publish_collection(self, collection_id):
        body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
        res = self.session.post(
            f"{self.api}/dp/v1/collections/{collection_id}/publish", headers=self.headers, data=json.dumps(body)
        )
        res.raise_for_status()


if __name__ == "__main__":
    smoke_test_init = SmokeTestsInitializer()
    # check whether we need to create collections
    test_collection_count = smoke_test_init.get_collection_count()
    visium_collection_count = smoke_test_init.get_collection_count(VISIUM_ACCT_CONTACT_NAME, 1)
    atac_seq_collection_count = smoke_test_init.get_collection_count(ATAC_SEQ_ACCT_CONTACT_NAME, 1)

    threads = []

    if test_collection_count < NUM_TEST_COLLECTIONS:
        num_to_create = NUM_TEST_COLLECTIONS - test_collection_count

        for _ in range(num_to_create):
            thread = threading.Thread(target=smoke_test_init.create_and_publish_collection)
            threads.append(thread)

    if visium_collection_count < 1:
        thread = threading.Thread(
            target=smoke_test_init.create_and_publish_collection,
            args=(VISIUM_ACCT_CONTACT_NAME, "Visium Test Collection", VISIUM_DATASET_MANIFEST),
        )
        threads.append(thread)

    if atac_seq_collection_count < 1:
        thread = threading.Thread(
            target=smoke_test_init.create_and_publish_collection,
            args=(ATAC_SEQ_ACCT_CONTACT_NAME, "ATAC Seq Test Collection", ATAC_SEQ_MANIFEST),
        )
        threads.append(thread)

    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
