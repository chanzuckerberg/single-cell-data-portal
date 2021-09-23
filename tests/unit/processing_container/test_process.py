import pathlib
from unittest.mock import patch

import requests

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
)
from backend.corpora.dataset_processing import process
from backend.corpora.dataset_processing.process import make_cxg, make_seurat, make_loom
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDatasetProcessing(CorporaTestCaseUsingMockAWS):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.real_h5ad_filename = pathlib.Path("fixtures/example_valid.h5ad").absolute()

    @staticmethod
    def get_presigned_url(dataset_id, asset_id):
        response = requests.post(
            f"https://api.cellxgene.staging.single-cell.czi.technology/dp/v1/datasets/" f"{dataset_id}/asset/{asset_id}"
        )
        response.raise_for_status()
        return response.json()["presigned_url"]

    @staticmethod
    def download(url, local_filename):
        with requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_filename, "wb") as fp:
                for chunk in resp.iter_content(chunk_size=None):
                    fp.write(chunk)

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    def test_make_cxg(self):
        make_cxg(str(self.real_h5ad_filename))

    def test_make_seurat(self):
        make_seurat(str(self.real_h5ad_filename))

    def test_make_loom(self):
        make_loom(str(self.real_h5ad_filename))

    @patch("backend.corpora.dataset_processing.process.download_from_dropbox_url")
    def test_main(self, mock_download_from_dropbox):
        """
        Tests full pipeline for processing an uploaded H5AD file, including
        file I/O and database updates, but excluding for the Dropbox download
        functionality.
        """
        mock_download_from_dropbox.return_value = self.real_h5ad_filename

        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )

        process.process(dataset.id, "https://www.dropbox.com/IGNORED", self.corpora_config.bucket_name, self.corpora_config.bucket_name)
