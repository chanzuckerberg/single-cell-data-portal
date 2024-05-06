import logging
import os
from typing import Union

from backend.cellguide.common.config import CellGuideConfig

logger = logging.getLogger(__name__)


def get_bucket_path() -> Union[str, None]:
    """
    This function generates a bucket path based on the deployment stage and remote development prefix.
    If the deployment stage is 'rdev' and a remote development prefix is provided, the bucket path is prefixed
    with 'env-rdev-cellguide' and the remote development prefix. Otherwise, the bucket path is the same as the
    bucket name.

    Returns
    -------
    str or None
        The generated bucket path if the deployment stage is set and valid, None otherwise.
    """

    deployment_stage = os.getenv("DEPLOYMENT_STAGE")
    remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX")

    if not deployment_stage:
        logger.warning("Not uploading the pipeline output to S3 as DEPLOYMENT_STAGE is not set.")
        return

    if deployment_stage == "rdev" and remote_dev_prefix is None:
        logger.warning(
            "Not uploading the pipeline output to S3 as REMOTE_DEV_PREFIX is not set when DEPLOYMENT_STAGE is rdev."
        )
        return
    elif deployment_stage == "rdev":
        bucket = CellGuideConfig().bucket
        bucket_path = f"s3://{bucket}/env-rdev-cellguide{remote_dev_prefix}/"
    elif deployment_stage in ["dev", "staging", "prod"]:
        bucket = CellGuideConfig().bucket
        bucket_path = f"s3://{bucket}/"
    else:
        logger.warning(
            f"Invalid DEPLOYMENT_STAGE value: {deployment_stage}. Please set DEPLOYMENT_STAGE to one of rdev, dev, staging, or prod"
        )
        return

    return bucket_path


def get_object_key(*, object: str) -> str:
    """
    This function generates an object key for S3 based on the deployment stage and remote development prefix.
    If the deployment stage is 'rdev' and a remote development prefix is provided, the object key is prefixed
    with 'env-rdev-cellguide' and the remote development prefix. Otherwise, the object key is the same as the
    input object.

    Parameters
    ----------
    object : str
        The input object for which the object key is to be generated.

    Returns
    -------
    str
        The generated object key.
    """

    deployment_stage = os.getenv("DEPLOYMENT_STAGE")
    remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX")

    object_key = object
    if deployment_stage == "rdev" and remote_dev_prefix is not None:
        CellGuideConfig().bucket
        object_key = f"env-rdev-cellguide{remote_dev_prefix}/{object}"

    return object_key
