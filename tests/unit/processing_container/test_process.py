import pathlib
import requests
import shutil
import tempfile

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
        cls.tmp_dir = tempfile.mkdtemp()
        cls.real_h5ad_filename = pathlib.Path(cls.tmp_dir, "real.h5ad")
        # TODO: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/corpora-data-portal/1354
        # Reupload this dataset once 2.0.0 migration is complete.
        cls.presigned_url = cls.get_presigned_url(
            "79dfc1fe-cc7a-454a-ba26-4c72f0876436", "c1f00b2e-fe42-4dd0-a035-7b7f3cefd055"
        )
        cls.download(cls.presigned_url, cls.real_h5ad_filename)

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
        shutil.rmtree(cls.tmp_dir)

    def test_make_cxg(self):
        make_cxg(str(self.real_h5ad_filename))

    def test_make_seurat(self):
        make_seurat(str(self.real_h5ad_filename))

    def test_make_loom(self):
        make_loom(str(self.real_h5ad_filename))

    def test_main(self):
        url = self.presigned_url
        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )
        self.bucket
        print("URL")
        print(url)
        process.process(dataset.id, url, self.corpora_config.bucket_name, self.corpora_config.bucket_name)
