import os
import shutil
from unittest.mock import patch

import requests

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
)
from backend.corpora.dataset_processing import process
from backend.corpora.dataset_processing.process import make_cxg, make_seurat, make_loom
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDatasetProcessing(CorporaTestCaseUsingMockAWS):
    @staticmethod
    def fixture_file_path(relative_filename):
        return os.path.abspath(os.path.join(os.path.dirname(__file__), relative_filename))

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.h5ad_with_labels = cls.fixture_file_path("fixtures/2_0_0_with_labels.h5ad")
        cls.h5ad_raw = cls.fixture_file_path("fixtures/2_0_0_raw_valid.h5ad")

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
        cls.clean_generated_files()

    @classmethod
    def clean_generated_files(cls):
        if os.path.exists(d := cls.fixture_file_path("local.cxg")):
            shutil.rmtree(d)
        if os.path.exists(f := cls.fixture_file_path("local.h5ad")):
            os.remove(f)
        if os.path.exists(d := cls.fixture_file_path("fixtures/2_0_0_with_labels.cxg")):
            shutil.rmtree(d)
        if os.path.exists(f := cls.fixture_file_path("fixtures/2_0_0_with_labels.rds")):
            os.remove(f)
        if os.path.exists(f := cls.fixture_file_path("fixtures/2_0_0_with_labels.loom")):
            os.remove(f)

    def test_make_cxg(self):
        make_cxg(str(self.h5ad_with_labels))

    def test_make_seurat(self):
        make_seurat(str(self.h5ad_with_labels))

    def test_make_loom(self):
        make_loom(str(self.h5ad_with_labels))

    @patch("backend.corpora.dataset_processing.process.download_from_dropbox_url")
    def test_main(self, mock_download_from_dropbox):
        """
        Tests full pipeline for processing an uploaded H5AD file, including
        file I/O and database updates, but excluding the Dropbox download
        functionality.
        """
        mock_download_from_dropbox.return_value = self.h5ad_raw

        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )

        process.process(dataset.id, "https://www.dropbox.com/IGNORED", self.corpora_config.bucket_name, self.corpora_config.bucket_name)

        # TODO: add assertions:
        # 1. H5AD has annotation labels added and uploaded to S3
        # 2. cxg, rds, loom uploaded to s3
        # 3. databases metadata updated and showing successful status
        # If we do the above, we can eliminate individual unit tests for #2, and the 2_0_0_with_labels.h5ad test fixture file
