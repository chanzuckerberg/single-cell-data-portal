#!/usr/bin/env python
import json
import os
import sys
import threading

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL
from tests.functional.backend.utils import (
    get_auth_token,
    make_cookie,
    make_proxy_auth_token,
    make_session,
    upload_and_wait,
)

# Amount to reduce chance of collision where multiple test instances select the same collection to test against
NUM_TEST_DATASETS = 3
NUM_TEST_COLLECTIONS = 10
TEST_ACCT_CONTACT_NAME = "Smoke Test User"


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
        self.headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

    def get_collection_count(self):
        res = self.session.get(f"{self.api}/curation/v1/collections?visiblity=PUBLIC", headers=self.headers)
        res.raise_for_status()
        data = json.loads(res.content)
        num_collections = 0
        for collection in data:
            if collection["contact_name"] == TEST_ACCT_CONTACT_NAME:
                num_collections += 1
            if num_collections == NUM_TEST_COLLECTIONS:
                return num_collections
        return num_collections

    def create_and_publish_collection(self, dropbox_url):
        collection_id = self.create_collection()
        for _ in range(NUM_TEST_DATASETS):
            upload_and_wait(self.session, self.api, self.curator_cookie, collection_id, dropbox_url)
        self.publish_collection(collection_id)
        print(f"created and published collection {collection_id}")

    def create_collection(self):
        data = {
            "contact_email": "example@gmail.com",
            "contact_name": TEST_ACCT_CONTACT_NAME,
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "http://protocol.com"}],
            "name": "test collection",
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
    collection_count = smoke_test_init.get_collection_count()
    if collection_count >= NUM_TEST_COLLECTIONS:
        sys.exit(0)
    dataset_dropbox_url = (
        "https://www.dropbox.com/scl/fi/y50umqlcrbz21a6jgu99z/"
        "5_0_0_example_valid.h5ad?rlkey=s7p6ybyx082hswix26hbl11pm&dl=0"
    )
    num_to_create = NUM_TEST_COLLECTIONS - collection_count
    threads = []
    for _ in range(num_to_create):
        thread = threading.Thread(target=smoke_test_init.create_and_publish_collection, args=(dataset_dropbox_url,))
        threads.append(thread)
        thread.start()
    for thread in threads:
        thread.join()
