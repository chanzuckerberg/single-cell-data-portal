import unittest
from unittest.mock import patch

from backend.common.feature_flag import FeatureFlagService, FeatureFlagValues


class TestFeatureFlag(unittest.TestCase):
    def setUp(self):
        super().setUp()

        def mock_config_fn(name):
            if name == FeatureFlagValues.SCHEMA_4:
                return "True"

        self.mock_config = patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
        self.mock_config.start()

    def tearDown(self):
        self.mock_config.stop()

    def test_feature_flag(self):
        self.assertTrue(FeatureFlagService.is_enabled(FeatureFlagValues.SCHEMA_4))
