from typing import Literal

from backend.common.corpora_config import CorporaConfig

"""
This is a lightweight wrapper around feature flags built with CorporaConfig. The intention is
to provide a generic interface for feature flags, so that if we ever want to swap out the
underlying implementation, we can do so without updating all client calls.

To add a new feature flag, add it to the FeatureFlag and FeatureFlagValues string literal above.

To use a feature flag:
```
if FeatureFlagService.is_enabled(FeatureFlagValues.<feature_flag_value_name>):
  <logic that only applies if feature flag is enabled>
```

To mock a feature flag in a test:
```
self.mock_config = CorporaConfig()
self.mock_config.set(dict(<feature_flag_value_name>="True"))
```
"""

FeatureFlag = Literal["citation_update_feature_flag"]


class FeatureFlagValues:
    CITATION_UPDATE = "citation_update_feature_flag"


class FeatureFlagService:
    @staticmethod
    def is_enabled(feature_flag: FeatureFlag) -> bool:
        flag_value = getattr(CorporaConfig(), feature_flag, "").lower()
        return flag_value == "true"
