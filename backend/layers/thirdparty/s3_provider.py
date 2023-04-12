import os
import subprocess
from typing import Tuple
from urllib.parse import urlparse

import boto3

from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface


class S3Provider(S3ProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("s3")

    def _parse_s3_uri(self, s3_uri: str) -> Tuple[str, str]:
        parsed_url = urlparse(s3_uri)
        return parsed_url.netloc, parsed_url.path[1:]

    def get_file_size(self, path: str) -> int:
        """
        Returns the file size of an S3 object located at `path`
        """
        bucket, key = self._parse_s3_uri(path)
        try:
            response = self.client.head_object(Bucket=bucket, Key=key)
        except Exception:
            return None
        return response["ContentLength"]

    def generate_presigned_url(self, path: str, expiration: int = 604800) -> str:
        """
        Generates a presigned url that can be used to download the S3 object located at `path`
        """
        bucket, key = self._parse_s3_uri(path)
        return self.client.generate_presigned_url(
            "get_object", Params={"Bucket": bucket, "Key": key}, ExpiresIn=expiration
        )

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        """
        Uploads the local `src_file` to an S3 object in the `bucket_name` bucket with object key `dst_file`
        """

        self.client.upload_file(
            src_file,
            bucket_name,
            dst_file,
            ExtraArgs=extra_args,
        )

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        """
        Downloads an S3 file located at s3://bucket_name/object_key to `local_filename`
        """
        self.client.download_file(bucket_name, object_key, local_filename)

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
