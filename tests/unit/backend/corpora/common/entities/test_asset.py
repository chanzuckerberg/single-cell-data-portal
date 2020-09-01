import unittest

from backend.corpora.common.entities.asset import Asset
from tests.unit.backend.corpora import CorporaTestCaseUsingMockAWS


class TestDataset(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_dataset_artifact_id"

    def test__get__ok(self):
        asset = Asset.get(self.uuid)
        self.assertEqual(self.uuid, asset.id)
        self.assertEqual(CorporaTestCaseUsingMockAWS.CORPORA_TEST_CONFIG["bucket_name"], asset.bucket_name)
        self.assertEqual("test_s3_uri.h5ad", asset.file_prefix)
