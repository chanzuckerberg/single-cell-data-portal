import json
import unittest
from unittest.mock import patch

from backend.common.corpora_orm import CollectionVisibility, ProcessingStatus
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatus, DatasetStatusKey
from tests.unit.backend.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.layers.common.base_api_test import DatasetStatusUpdate


@patch(
    "backend.common.utils.dl_sources.url.DropBoxURL.file_info",
    return_value={"size": 1, "name": "file.h5ad"},
)
@patch("backend.common.upload.start_upload_sfn")
class TestPutLink(BaseAuthAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        cls.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    def test__from_link__no_auth(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 401 Unauthorized if the user is not authenticated
        """
        dataset = self.generate_dataset(statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)])
        body = {"link": self.good_link}
        headers = None
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )

        self.assertEqual(401, response.status_code)

    def test__from_link__Not_Public(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 403 Unauthorized
        if trying to upload to a published collection
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
            publish=True,
        )
        body = {"link": self.good_link}
        headers = self.make_owner_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )

        self.assertEqual(403, response.status_code)

    def test__from_link__Not_Owner(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 403 Unauthorized
        if the authenticated user is not the owner of a collection
        """

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
        )
        body = {"link": self.dummy_link}
        headers = self.make_not_owner_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )

        self.assertEqual(403, response.status_code)

    def test__new_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should succeed if a valid link is uploaded by the owner of the collection
        """

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
        )
        body = {"link": self.good_link}
        headers = self.make_owner_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )
        self.assertEqual(202, response.status_code)

    def test__new_from_link__Super_Curator(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should succeed if a valid link is uploaded by a super curator
        """

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
        )
        body = {"link": self.good_link}
        headers = self.make_super_curator_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )
        self.assertEqual(202, response.status_code)

    def test__existing_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id on an existing dataset_id 
        should succeed if a valid link is uploaded by the owner of the collection
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )
        body = {"link": self.good_link}
        headers = self.make_owner_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )
        self.assertEqual(202, response.status_code)

    def test__existing_from_link__Super_Curator(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id on an existing dataset_id 
        should succeed if a valid link is uploaded by a super curator
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )
        body = {"link": self.good_link}
        headers = self.make_super_curator_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}", json=body, headers=headers
        )
        self.assertEqual(202, response.status_code)


if __name__ == "__main__":
    unittest.main()
