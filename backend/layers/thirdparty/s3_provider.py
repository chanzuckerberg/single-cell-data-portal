import logging
import os
import re
import subprocess
from typing import Iterable, List, Tuple
from urllib.parse import urlparse

import boto3

from backend.layers.common.regex import ID_REGEX
from backend.layers.thirdparty.s3_exceptions import IllegalS3RecursiveDelete, S3DeleteException
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface

AWS_S3_MAX_ITEMS_PER_BATCH = 1000

logger = logging.getLogger(__name__)


def _parse_size_string(size_str: str) -> int:
    """
    Parse a human-readable size string (e.g., "50GB", "100MB", "1TB") into bytes.

    Supports: TB, GB, MB, KB, B (case-insensitive)
    Uses binary units (1024-based).

    Args:
        size_str: Size string in format like "50GB", "100MB", etc.

    Returns:
        Size in bytes as integer

    Raises:
        ValueError: If the format is invalid
    """
    if not size_str:
        raise ValueError("Size string cannot be empty")

    # Match number and unit (case-insensitive)
    # Pattern matches: number (int or float) + optional whitespace + unit (TB, GB, MB, KB, or B)
    match = re.match(r"^(\d+(?:\.\d+)?)\s*((?:TB|GB|MB|KB|B))$", size_str.strip(), re.IGNORECASE)
    if not match:
        raise ValueError(f"Invalid size format: {size_str}. Expected format: '50GB', '100MB', etc.")

    value = float(match.group(1))
    unit = match.group(2).upper()

    # Convert to bytes based on unit (binary units: 1024-based)
    multipliers = {
        "TB": 1024**3,  # 1 TB = 1024³ bytes
        "GB": 1024**2,  # 1 GB = 1024² bytes
        "MB": 1024**1,  # 1 MB = 1024 bytes
        "KB": 1024,  # 1 KB = 1024 bytes
        "B": 1,  # 1 B = 1 byte
    }

    if unit not in multipliers:
        raise ValueError(f"Unsupported unit: {unit}. Supported units: TB, GB, MB, KB, B")

    return int(value * multipliers[unit])


def _get_cloudfront_max_cache_size() -> int:
    """
    Get the CloudFront max cache size from environment variable or use default.

    Reads CLOUDFRONT_MAX_CACHE_SIZE environment variable and parses it.
    If not set or invalid, defaults to 50GB.

    Returns:
        Max cache size in bytes
    """
    default_size_bytes = 50 * 1024 * 1024 * 1024  # 50GB

    env_value = os.getenv("CLOUDFRONT_MAX_CACHE_SIZE")
    if not env_value:
        return default_size_bytes

    try:
        return _parse_size_string(env_value)
    except ValueError as e:
        logger.warning(
            {
                "message": "Invalid CLOUDFRONT_MAX_CACHE_SIZE format, using default 50GB",
                "env_value": env_value,
                "error": str(e),
            }
        )
        return default_size_bytes


CLOUDFRONT_MAX_CACHE_SIZE = _get_cloudfront_max_cache_size()


class S3Provider(S3ProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("s3")

    @staticmethod
    def parse_s3_uri(s3_uri: str) -> Tuple[str, str]:
        parsed_url = urlparse(s3_uri)
        return parsed_url.netloc, parsed_url.path[1:]

    def uri_exists(self, s3_uri: str) -> bool:
        bucket, key = self.parse_s3_uri(s3_uri)
        try:
            logger.info(
                {
                    "message": "Running HEAD request",
                    "s3_uri": s3_uri,
                    "bucket": bucket,
                    "key": key,
                }
            )
            self.client.head_object(Bucket=bucket, Key=key)
            logger.info(
                {
                    "message": "HEAD request succeeded!",
                }
            )
            return True
        except Exception as e:
            logger.info(
                {
                    "message": "HEAD request failed",
                    "error": str(e),
                }
            )
            return False

    def get_file_size(self, path: str) -> int:
        """
        Returns the file size of an S3 object located at `path`
        """
        bucket, key = self.parse_s3_uri(path)
        try:
            response = self.client.head_object(Bucket=bucket, Key=key)
        except Exception:
            return None
        return response["ContentLength"]

    def generate_presigned_url(self, path: str, expiration: int = 604800) -> str:
        """
        Generates a presigned url that can be used to download the S3 object located at `path`
        """
        bucket, key = self.parse_s3_uri(path)
        return self.client.generate_presigned_url(
            "get_object", Params={"Bucket": bucket, "Key": key}, ExpiresIn=expiration
        )

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        """
        Uploads the local `src_file` to an S3 object in the `bucket_name` bucket with object key `dst_file`
        Sets Cache-Control: no-cache for files larger than CloudFront max cache size to prevent CloudFront caching issues
        """
        # Check file size and set Cache-Control for large files (exceeding CloudFront max cache size)
        file_size = os.path.getsize(src_file)

        # Create a copy of extra_args to avoid modifying the original
        upload_args = extra_args.copy() if extra_args else {}

        if file_size > CLOUDFRONT_MAX_CACHE_SIZE:
            # Set Cache-Control to no-cache to bypass CloudFront caching for files larger than max cache size
            upload_args["CacheControl"] = "no-cache, no-store, must-revalidate"
            logger.info(
                {
                    "message": "Large file detected - setting Cache-Control: no-cache",
                    "file_size_gb": file_size / (1024 * 1024 * 1024),
                    "src_file": src_file,
                    "dst_file": dst_file,
                }
            )

        logger.info(
            {
                "message": "Uploading file",
                "bucket_name": bucket_name,
                "dst_file": dst_file,
                "src_file": src_file,
                "file_size_bytes": file_size,
                "extra_args": upload_args,
            }
        )
        self.client.upload_file(
            src_file,
            bucket_name,
            dst_file,
            ExtraArgs=upload_args,
        )

    def delete_files(self, bucket_name: str, object_keys: List[str]) -> None:
        """
        Deletes the objects `object_keys` from bucket `bucket_name`
        """
        for i in range(0, len(object_keys), AWS_S3_MAX_ITEMS_PER_BATCH):
            key_batch = object_keys[i : i + AWS_S3_MAX_ITEMS_PER_BATCH]
            resp = self.client.delete_objects(
                Bucket=bucket_name,
                Delete={"Objects": [{"Key": key} for key in key_batch]},
            )
            if deleted := resp.get("Deleted"):
                logger.info({"deleted": deleted})
            if errors := resp.get("Errors"):
                logger.error({"errors": errors, "bucket_name": bucket_name})
                raise S3DeleteException(errors)

    def delete_prefix(self, bucket_name: str, prefix: str) -> None:
        if not re.search(ID_REGEX, prefix):
            raise IllegalS3RecursiveDelete("Cannot recursively delete without a valid UUID prefix")
        object_keys = list(self.list_directory(bucket_name, prefix))
        self.delete_files(bucket_name, object_keys)

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        """
        Downloads an S3 file located at s3://bucket_name/object_key to `local_filename`
        """
        logger.info({"message": "Downloading file", "bucket_name": bucket_name, "object_key": object_key})
        self.client.download_file(bucket_name, object_key, local_filename)

    def restore_object(self, bucket_name: str, object_key: str) -> None:
        response = self.client.list_object_versions(
            Bucket=bucket_name,
            Prefix=object_key,
        )
        if delete_markers := response.get("DeleteMarkers"):
            for marker in delete_markers:
                if not marker["IsLatest"]:
                    continue
                logger.info("restoring", marker["Key"])
                response = self.client.delete_object(
                    Bucket=bucket_name,
                    Key=marker["Key"],
                    VersionId=marker["VersionId"],
                )
                response.pop("ResponseMetadata")
                logger.info(response)

    def upload_directory(self, src_dir: str, s3_uri: str):
        """
        Uploads a whole local directory `src_dir` to `s3_uri`
        """
        command = ["aws"]
        if os.getenv("BOTO_ENDPOINT_URL"):
            command.append(f"--endpoint-url={os.getenv('BOTO_ENDPOINT_URL')}")

        logger.info(
            {
                "message": "Copying directory",
                "src_dir": src_dir,
                "dst_dir": s3_uri,
            }
        )

        command.extend(
            [
                "s3",
                "cp",
                src_dir,
                s3_uri,
                "--recursive",
                "--acl",
                "bucket-owner-full-control",
            ]
        )
        subprocess.run(
            command,
            check=True,
        )

    def list_directory(self, bucket_name: str, src_dir: str) -> Iterable[str]:
        paginator = self.client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket_name, Prefix=src_dir):
            yield from (content["Key"] for content in page.get("Contents", []))

    def put_object(self, bucket_name: str, object_key: str, body: str) -> None:
        """
        Uploads a json object to S3
        """
        self.client.put_object(
            Bucket=bucket_name,
            Key=object_key,
            Body=body,
        )

    def create_bucket(self, bucket_name: str, location: str) -> None:
        """
        Creates a bucket
        """
        self.client.create_bucket(Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": f"{location}"})
