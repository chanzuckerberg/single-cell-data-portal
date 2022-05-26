import json
import unittest
from unittest.mock import patch

from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestPutLink(BaseAuthAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        cls.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    def _test_new(self, collection_params: dict = None, headers: dict = None, body: dict = None):
        headers = headers if headers else {}
        collection_params = collection_params if collection_params else {}
        body = body if body else ""
        headers["Content-Type"] = "application/json"
        collection = self.generate_collection(self.session, **collection_params)
        response = self.app.put(
            f"/curation/v1/collections/{collection.id}/datasets/upload-link", data=json.dumps(body), headers=headers
        )
        return response

    def test__from_link__no_auth(self):
        response = self._test_new()
        self.assertEqual(401, response.status_code)

    def test__from_link__Not_Public(self):
        response = self._test_new(
            dict(visibility=CollectionVisibility.PUBLIC.name),
            self.get_auth_headers(),
            body={"curator_tag": "test", "link": self.dummy_link},
        )
        self.assertEqual(403, response.status_code)

    def test__from_link__Not_Owner(self):
        response = self._test_new(
            dict(owner="someone else"), self.get_auth_headers(), body={"curator_tag": "test", "link": self.dummy_link}
        )
        self.assertEqual(403, response.status_code)

    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__new_from_link__OK(self, start_upload_sfn):
        response = self._test_new({}, self.get_auth_headers(), body={"curator_tag": "test", "link": self.good_link})
        self.assertEqual(202, response.status_code)

    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__new_from_link__Super_Curator(self, start_upload_sfn):
        headers = self.make_super_curator_header()
        response = self._test_new({}, headers, body={"curator_tag": "test", "link": self.good_link})
        self.assertEqual(202, response.status_code)

    def _test_existing(self, headers: dict = None, use_curator_tag=False):
        curator_tag = "test.h5ad"
        headers = headers if headers else {}
        headers["Content-Type"] = "application/json"
        collection = self.generate_collection(self.session)
        processing_status = dict(processing_status=ProcessingStatus.SUCCESS)
        dataset = self.generate_dataset(
            self.session, collection_id=collection.id, curator_tag=curator_tag, processing_status=processing_status
        )
        body = {"id": dataset.id, "link": self.good_link}
        if use_curator_tag:
            body["curator_tag"] = curator_tag
        response = self.app.put(
            f"/curation/v1/collections/{collection.id}/datasets/upload-link", data=json.dumps(body), headers=headers
        )
        return response

    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__existing_from_link__OK(self, start_upload_sfn):
        with self.subTest("dataset_id"):
            response = self._test_existing(self.get_auth_headers())
            self.assertEqual(202, response.status_code)
        with self.subTest("curator_tag"):
            response = self._test_existing(self.get_auth_headers(), use_curator_tag=True)
            self.assertEqual(202, response.status_code)

    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__existing_from_link__Super_Curator(self, start_upload_sfn):
        headers = self.make_super_curator_header()
        with self.subTest("dataset_id"):
            response = self._test_existing(headers, use_curator_tag=True)
            self.assertEqual(202, response.status_code)
        with self.subTest("curator_tag"):
            response = self._test_existing(headers, use_curator_tag=True)
            self.assertEqual(202, response.status_code)

    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__curator_tag_ignored_when_dataset_id_is_present__OK(self, start_upload_sfn):
        curator_tag = "test.h5ad"
        headers = self.get_auth_headers()
        headers["Content-Type"] = "application/json"
        collection = self.generate_collection(self.session)
        processing_status = dict(processing_status=ProcessingStatus.SUCCESS)
        dataset = self.generate_dataset(
            self.session, collection_id=collection.id, curator_tag=curator_tag, processing_status=processing_status
        )

        body = {"id": dataset.id, "link": self.good_link, "curator_tag": "different_tag"}
        response = self.app.put(
            f"/curation/v1/collections/{collection.id}/datasets/upload-link", data=json.dumps(body), headers=headers
        )
        self.assertEqual(202, response.status_code)


if __name__ == "__main__":
    unittest.main()
