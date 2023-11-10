from enum import Enum

from backend.common.corpora_config import CorporaConfig

"""
This is a lightweight wrapper around feature flags built with CorporaConfig. The intention is
to provide a generic interface for feature flags, so that if we ever want to swap out the
underlying implementation, we can do so without updating all client calls.

To add a new feature flag, add it to the FeatureFlag and FeatureFlag string literal above.

To use a feature flag:
```
if FeatureFlagService.is_enabled(FeatureFlag.SCHEMA_4):
  <logic that only applies if schema 4 is enabled>
```

To mock a feature flag in a test:
```
def mock_config_fn(name):
    if name == FeatureFlag.SCHEMA_4:
        return "True"

self.mock_config = patch.object(FeatureFlagService, "is_enabled", side_effect=mock_config_fn)
```
"""


class FeatureFlag(Enum):
    SCHEMA_4 = "schema_4_feature_flag"


class FeatureFlagService:
    @staticmethod
    def is_enabled(feature_flag: FeatureFlag) -> bool:
        flag_value = getattr(CorporaConfig(), feature_flag.value, "").lower()
        return flag_value == "true"
