import json
import unittest
from unittest.mock import patch

from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


@patch(
    "backend.corpora.common.utils.dl_sources.url.DropBoxURL.file_info",
    return_value={"size": 1, "name": "file.h5ad"},
)
@patch("backend.corpora.common.upload.start_upload_sfn")
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

    def test__from_link__no_auth(self, *mocks):
        response = self._test_new()
        self.assertEqual(401, response.status_code)

    def test__from_link__Not_Public(self, *mocks):
        response = self._test_new(
            dict(visibility=CollectionVisibility.PUBLIC.name),
            self.make_owner_header(),
            body={"link": self.dummy_link},
        )
        self.assertEqual(403, response.status_code)

    def test__from_link__Not_Owner(self, *mocks):
        response = self._test_new(
            dict(owner="someone else"),
            self.make_not_owner_header(),
            body={"link": self.dummy_link},
        )
        self.assertEqual(403, response.status_code)

    def test__new_from_link__OK(self, *mocks):
        headers = self.make_owner_header()
        response = self._test_new({}, headers, body={"link": self.good_link})
        self.assertEqual(202, response.status_code)

    def test__new_from_link__Super_Curator(self, *mocks):
        headers = self.make_super_curator_header()
        response = self._test_new({}, headers, body={"link": self.good_link})
        self.assertEqual(202, response.status_code)

    def _test_existing(self, headers: dict = None, use_curator_tag=False):
        headers = headers if headers else {}
        headers["Content-Type"] = "application/json"
        collection = self.generate_collection(self.session)
        processing_status = dict(processing_status=ProcessingStatus.SUCCESS)
        dataset = self.generate_dataset(self.session, collection_id=collection.id, processing_status=processing_status)
        body = {"id": dataset.id, "link": self.good_link}
        response = self.app.put(
            f"/curation/v1/collections/{collection.id}/datasets/upload-link", data=json.dumps(body), headers=headers
        )
        return response

    def test__existing_from_link__OK(self, *mocks):
        with self.subTest("dataset_id"):
            response = self._test_existing(self.make_owner_header())
            self.assertEqual(202, response.status_code)

    def test__existing_from_link__Super_Curator(self, *mocks):
        headers = self.make_super_curator_header()
        with self.subTest("dataset_id"):
            response = self._test_existing(headers, use_curator_tag=True)
            self.assertEqual(202, response.status_code)


if __name__ == "__main__":
    unittest.main()
