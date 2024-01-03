import unittest
from unittest.mock import Mock, patch

from backend.common.corpora_config import CorporaConfig
from backend.common.feature_flag import FeatureFlagService


class TestFeatureFlag(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.mock_config = CorporaConfig()
        self.mock_config.set({"test_feature_flag": "True"})

    def tearDown(self):
        self.mock_config.reset()

    @patch("backend.common.feature_flag.FeatureFlag", Mock())
    def test_feature_flag(self):
        self.assertTrue(FeatureFlagService.is_enabled("test_feature_flag"))
