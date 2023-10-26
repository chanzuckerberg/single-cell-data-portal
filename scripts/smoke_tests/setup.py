#!/usr/bin/env python
import json
import sys
import threading

from tests.functional.backend.common import BaseFunctionalTestCase

# Amount to reduce chance of collision where multiple test instances select the same collection to test against
NUM_TEST_COLLECTIONS = 10
TEST_ACCT_CONTACT_NAME = "Smoke Test User"


class SmokeTestsInitializer(BaseFunctionalTestCase):
    def __init__(self):
        super().setUpClass(smoke_tests=True)
        super().setUp()
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
        self.upload_and_wait(collection_id, dropbox_url, cleanup=False)
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
        body = {"data_submission_policy_version": "1.0"}
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
    if smoke_test_init.is_using_schema_4:
        dataset_dropbox_url = (
            "https://www.dropbox.com/scl/fi/d99hpw3p2cxtmi7v4kyv5/"
            "4_0_0_test_dataset.h5ad?rlkey=i5ownt8g1mropbu41r7fa0i06&dl=0"
        )
    else:
        dataset_dropbox_url = "https://www.dropbox.com/s/m1ur46nleit8l3w/3_0_0_valid.h5ad?dl=0"
    num_to_create = NUM_TEST_COLLECTIONS - collection_count
    threads = []
    for _ in range(num_to_create):
        thread = threading.Thread(target=smoke_test_init.create_and_publish_collection, args=(dataset_dropbox_url,))
        threads.append(thread)
        thread.start()
    for thread in threads:
        thread.join()
