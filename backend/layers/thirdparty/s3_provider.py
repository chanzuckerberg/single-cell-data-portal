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


class S3Provider(S3ProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("s3")

    @staticmethod
    def parse_s3_uri(s3_uri: str) -> Tuple[str, str]:
        parsed_url = urlparse(s3_uri)
        return parsed_url.netloc, parsed_url.path[1:]

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
        """
        logger.info(
            {
                "message": "Uploading file",
                "bucket_name": bucket_name,
                "dst_file": dst_file,
                "src_file": src_file,
                "extra_args": extra_args,
            }
        )
        self.client.upload_file(
            src_file,
            bucket_name,
            dst_file,
            ExtraArgs=extra_args,
        )

    def delete_files(self, bucket_name: str, object_keys: List[str]) -> None:
        """
        Deletes the objects `object_keys` from bucket `bucket_name`
        """
        for i in range(0, len(object_keys), AWS_S3_MAX_ITEMS_PER_BATCH):
            key_batch = object_keys[i:AWS_S3_MAX_ITEMS_PER_BATCH]
            resp = self.client.delete_objects(
                Bucket=bucket_name,
                Delete={"Objects": [{"Key": key} for key in key_batch]},
            )
            if deleted := resp.get("Deleted"):
                logger.info(f"Deleted: {deleted}")
            if errors := resp.get("Errors"):
                logger.info(f"Errors: {errors}")
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
