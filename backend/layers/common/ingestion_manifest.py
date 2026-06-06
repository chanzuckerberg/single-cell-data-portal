import json
import re
from typing import Optional, Union

import pydantic
from pydantic import AnyUrl, BaseModel, HttpUrl

_PYDANTIC_V1 = int(pydantic.VERSION.split(".")[0]) < 2


class S3Url(AnyUrl):
    """Pydantic Model for S3 URLs

    Copied from https://gist.github.com/rajivnarayan/c38f01b89de852b3e7d459cfde067f3f
    # TODO consolidate with backend/common/utils/dl_sources/uri.py
    """

    allowed_schemes = {"s3"}
    pattern = re.compile(
        r"^s3://"
        r"(?=[a-z0-9])"  # Bucket name must start with a letter or digit
        r"(?!(^xn--|sthree-|sthree-configurator|.+-s3alias$))"  # Bucket name must not start with xn--, sthree-, sthree-configurator or end with -s3alias
        r"(?!.*\.\.)"  # Bucket name must not contain two adjacent periods
        r"[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]"  # Bucket naming constraints
        r"(?<!\.-$)"  # Bucket name must not end with a period followed by a hyphen
        r"(?<!\.$)"  # Bucket name must not end with a period
        r"(?<!-$)"  # Bucket name must not end with a hyphen
        r"(/([a-zA-Z0-9._-]+/?)*)?$"  # key naming constraints
    )


class IngestionManifest(BaseModel):
    """
    # Deserialize JSON to Pydantic model
    data_obj = IngestionManifest.model_validate_json(json_data)

    # Convert Pydantic object to JSON
    json_output = data_obj.model_dump_json(indent=2)
    """

    anndata: Union[HttpUrl, S3Url]
    atac_fragment: Optional[Union[HttpUrl, S3Url]] = None  # Optional field
    is_pre_analysis: bool = False


if _PYDANTIC_V1:
    # Add pydantic v2-compatible methods when running under pydantic v1.
    # This allows the codebase to use the pydantic v2 API uniformly.

    def _model_dump(self, exclude_none: bool = False, **kwargs):
        raw = self.dict(**kwargs)
        if exclude_none:
            raw = {k: v for k, v in raw.items() if v is not None}
        return raw

    def _model_dump_json(self, **kwargs) -> str:
        # Serialize URL types to strings, keep None/bool as-is for valid JSON output
        d = {}
        for k, v in self.dict().items():
            if isinstance(v, bool) or v is None:
                d[k] = v
            else:
                d[k] = str(v)
        return json.dumps(d, separators=(",", ":"))

    @classmethod
    def _model_validate_json(cls, json_str: str) -> "IngestionManifest":
        return cls.parse_raw(json_str)

    IngestionManifest.model_dump = _model_dump
    IngestionManifest.model_dump_json = _model_dump_json
    IngestionManifest.model_validate_json = classmethod(_model_validate_json.__func__)
