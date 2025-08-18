import contextlib
import logging
import re
from typing import Iterable, List, Set
from urllib.parse import urlparse

from backend.layers.common.regex import ID_REGEX
from backend.layers.thirdparty.s3_exceptions import IllegalS3RecursiveDelete
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface

logger = logging.getLogger(__name__)


class MockS3Provider(S3ProviderInterface):
    """
    Simple S3 mock that mostly checks if paths are correct
    """

    def __init__(self) -> None:
        self.mock_s3_fs: Set[str] = set()

    def parse_s3_uri(self, s3_uri: str):
        parsed_url = urlparse(s3_uri)
        return parsed_url.netloc, parsed_url.path[1:]

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        url = f"s3://{bucket_name}/{dst_file}"
        self.mock_s3_fs.add(url)

    def upload_directory(self, src_dir: str, s3_uri: str):
        self.mock_s3_fs.add(s3_uri)

    def delete_files(self, bucket_name: str, object_keys: List[str]):
        # Filter out empty or invalid keys
        valid_keys = [key.strip() for key in object_keys if key and key.strip()]

        if not valid_keys:
            logger.info({"message": "No valid keys to delete", "bucket_name": bucket_name})
            return

        # Safety check: prevent deletion of root-level objects without proper prefixes
        # Allow UUID-based filenames (36 chars) for public bucket compatibility
        dangerous_keys = []
        for key in valid_keys:
            if "/" not in key:
                # Root-level file - check if it's a UUID-based filename
                filename_without_ext = key.split(".")[0]
                if not re.match(ID_REGEX, filename_without_ext):
                    dangerous_keys.append(key)
            elif len(key.split("/")[0]) < 8:
                # Directory structure but prefix too short
                dangerous_keys.append(key)

        if dangerous_keys:
            logger.warning(
                {
                    "message": "Blocked deletion of potentially dangerous keys without proper directory structure or UUID filename",
                    "bucket_name": bucket_name,
                    "dangerous_keys": dangerous_keys[:10],
                }
            )
            raise IllegalS3RecursiveDelete(
                f"Cannot delete root-level or insufficiently prefixed objects: {dangerous_keys[:5]}"
            )

        # Safety check: warn about large deletion operations
        if len(valid_keys) > 10000:
            logger.warning(
                {
                    "message": "Large deletion operation detected - proceeding with caution",
                    "bucket_name": bucket_name,
                    "key_count": len(valid_keys),
                    "sample_keys": valid_keys[:5],
                }
            )

        # Log all deletion operations for audit trail
        logger.info(
            {
                "message": "Starting deletion operation",
                "bucket_name": bucket_name,
                "key_count": len(valid_keys),
                "sample_keys": valid_keys[:3],
            }
        )

        for key in valid_keys:
            with contextlib.suppress(KeyError):  # if key is not in bucket
                self.mock_s3_fs.remove(f"s3://{bucket_name}/{key}")

    def delete_prefix(self, bucket_name: str, prefix: str) -> None:
        object_keys = list(self.list_directory(bucket_name, prefix))
        self.delete_files(bucket_name, object_keys)

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        pass

    def restore_object(self, bucket_name: str, object_key: str) -> None:
        url = f"s3://{bucket_name}/{object_key}"
        self.mock_s3_fs.add(url)

    def file_exists(self, bucket_name: str, object_key: str):
        url = f"s3://{bucket_name}/{object_key}"
        return self.uri_exists(url)

    def list_directory(self, bucket_name: str, src_dir: str) -> Iterable[str]:
        object_prefix = f"s3://{bucket_name}/{src_dir}"
        for key in self.mock_s3_fs:
            if key.startswith(object_prefix):
                prefix_index = key.find(bucket_name) + len(bucket_name + "/")
                yield key[prefix_index:]

    def uri_exists(self, uri: str):
        return uri in self.mock_s3_fs

    def is_empty(self):
        return len(self.mock_s3_fs) == 0
