import os
import re

from backend.common.utils.regex import DATASET_ID_REGEX


def rds_citation_from_h5ad(citation: str) -> str:
    return re.sub(f"{DATASET_ID_REGEX}\.h5ad", r"\g<dataset_id>.rds", citation)


def get_s3_key_for_migration_record(
    schema_version: str, execution_id: str, dataset_id: str, old_version_id: str, new_version_id: str, file_name: str
) -> str:
    remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "").strip("/")
    if remote_dev_prefix:
        remote_dev_prefix = remote_dev_prefix + "/"
    return (
        f"{remote_dev_prefix}schema_migration/{schema_version}/{execution_id}/processed/dataset_{dataset_id}"
        f"/old_version_{old_version_id}_new_version_{new_version_id}/{file_name}"
    )
