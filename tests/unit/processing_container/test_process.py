import os
import shutil
from unittest.mock import patch

import requests

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
)
from backend.corpora.dataset_processing import process
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDatasetProcessing(CorporaTestCaseUsingMockAWS):
    """
    This test class is intended to exercise the upload pipeline processing, running within
    a "corpora-upload" (upload processing) Docker container. As the Docker container is configured
    with all necessary software dependencies, manually running this test outside of Docker
    (e.g. within an IDE when making changes to the test itself) will require you to:
    * set env vars: CORPORA_LOCAL_DEV=1;BOTO_ENDPOINT_URL=http://localhost:4566
    * Locally install latest cellxgene-schema (via pip)
    * Install all the same Python dependencies at Dockerfile.processing_image:23 and Dockerfile.processing_base:19
      (R software) to ensure Loom and Seurat artifacts can be generated. (These will only matter if we add test
      assertions to check these artifact's creation.)
    """

    @staticmethod
    def fixture_file_path(relative_filename):
        return os.path.abspath(os.path.join(os.path.dirname(__file__), relative_filename))

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
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
        if os.path.exists("local.cxg"):
            shutil.rmtree("local.cxg")
        for f in ["local.h5ad", "local.loom", "local.rds"]:
            if os.path.exists(f):
                os.remove(f)

    @patch("backend.corpora.dataset_processing.process.download_from_dropbox_url")
    def test_main(self, mock_download_from_dropbox):
        """
        Tests full pipeline for processing an uploaded H5AD file, including database updates
        generation and upload of all artifacts to S3 (localstack), but excluding the Dropbox download
        functionality.  Dropbox I/O is is mocked to prevent dependency on remote services (non-Dockerized).
        """
        mock_download_from_dropbox.return_value = self.h5ad_raw

        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )

        process.process(
            dataset.id,
            "https://www.dropbox.com/IGNORED",
            self.corpora_config.bucket_name,
            self.corpora_config.bucket_name,
        )

        # TODO: add assertions:
        # 1. H5AD has annotation labels added and uploaded to S3
        # 2. cxg, rds, loom uploaded to s3
        # 3. databases metadata updated and showing successful status
