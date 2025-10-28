import re
from typing import Optional, Union

from pydantic import AnyUrl, BaseModel, HttpUrl


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
