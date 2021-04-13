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
        presigned_url = cls.get_presigned_url(
            "a5da03d7-f44b-4bb8-81a4-bc80960d6223", "04be802c-ffc6-4a7f-a4c5-f2812f1a707f"
        )
        cls.download(presigned_url, cls.real_h5ad_filename)

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
        with self.assertRaises(RuntimeError):
            make_seurat(str(self.real_h5ad_filename))

    def test_make_loom(self):
        make_loom(str(self.real_h5ad_filename))

    def test_main(self):
        url = self.get_presigned_url("a5da03d7-f44b-4bb8-81a4-bc80960d6223", "04be802c-ffc6-4a7f-a4c5-f2812f1a707f")
        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )
        self.bucket
        process.process(dataset.id, url, self.corpora_config.bucket_name, self.corpora_config.bucket_name)
