import unittest

from backend.common.corpora_config import CorporaConfig
from backend.common.feature_flag import FeatureFlagService, FeatureFlagValues


class TestFeatureFlag(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.mock_config = CorporaConfig()
        self.mock_config.set({"schema_4_feature_flag": "True"})

    def tearDown(self):
        self.mock_config.reset()

    def test_feature_flag(self):
        self.assertTrue(FeatureFlagService.is_enabled(FeatureFlagValues.SCHEMA_4))
