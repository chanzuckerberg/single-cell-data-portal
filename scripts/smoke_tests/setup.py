#!/usr/bin/env python
import sys
import json
import time
import threading

import requests

from tests.functional.backend.common import BaseFunctionalTestCase

# Amount to reduce chance of collision where multiple test instances select the same collection to test against
NUM_TEST_COLLECTIONS = 10
TEST_ACCT_CONTACT_NAME = "Smoke Test User"


class SmokeTestsInitializer(BaseFunctionalTestCase):

    def __init__(self):
        super().setUpClass()
        super().setUp()
        self.headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

    def get_collection_count(self):
        res = self.session.get(f"{self.api}/curation/v1/collections?visiblity=PUBLIC", headers=self.headers)
        res.raise_for_status()
        data = json.loads(res.content)
        num_collections = 0
        for collection in data['collections']:
            if collection["contact_name"] == TEST_ACCT_CONTACT_NAME:
                num_collections += 1
            if num_collections == NUM_TEST_COLLECTIONS:
                return num_collections
        return num_collections

    def create_and_publish_collection(self, dropbox_url):
        collection_id = self.create_collection()
        self.upload_and_wait(collection_id, dropbox_url)
        self.publish_collection(collection_id)
        print(f"created and published collection {collection_id}")

    def create_collection(self):
        data = {
            "contact_email": "example@gmail.com",
            "contact_name": TEST_ACCT_CONTACT_NAME,
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "http://protocol.com"}],
            "name": "test collection",
        }

        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=self.headers)
        res.raise_for_status()
        data = json.loads(res.content)
        return data["collection_id"]

    def publish_collection(self, collection_id):
        body = {"data_submission_policy_version": "1.0"}
        res = self.session.post(
            f"{self.api}/dp/v1/collections/{collection_id}/publish", headers=self.headers, data=json.dumps(body)
        )
        res.raise_for_status()

    # override
    def upload_and_wait(self, collection_id, dropbox_url):
        body = {"url": dropbox_url}
        res = self.session.post(
            f"{self.api}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=self.headers
        )
        res.raise_for_status()
        dataset_id = json.loads(res.content)["dataset_id"]
        keep_trying = True
        timer = time.time()
        while keep_trying:
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=self.headers)
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
            time.sleep(20)
        return dataset_id


if __name__ == "__main__":
    smoke_test_init = SmokeTestsInitializer()
    # check whether we need to create collections
    collection_count = smoke_test_init.get_collection_count()
    if collection_count >= NUM_TEST_COLLECTIONS:
        print("Found sufficient published collections for testing, exiting")
        sys.exit(0)

    dataset_dropbox_url = "https://www.dropbox.com/s/qiclvn1slmap351/example_valid.h5ad?dl=0"
    num_to_create = NUM_TEST_COLLECTIONS - collection_count
    threads = []
    for i in range(num_to_create):
        thread = threading.Thread(target=smoke_test_init.create_and_publish_collection, args=(dataset_dropbox_url,))
        threads.append(thread)
        thread.start()
    for thread in threads:
        thread.join()