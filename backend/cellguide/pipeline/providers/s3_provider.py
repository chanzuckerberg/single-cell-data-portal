import logging
import os
import subprocess

import boto3
import botocore

AWS_S3_MAX_ITEMS_PER_BATCH = 1000

logger = logging.getLogger(__name__)


class S3Provider:
    def __init__(self) -> None:
        self.client = boto3.client("s3")

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

    def sync_directory(self, src_dir: str, s3_uri: str):
        """
        Uploads a whole local directory `src_dir` to `s3_uri`
        """
        command = ["aws"]
        if os.getenv("BOTO_ENDPOINT_URL"):
            command.append(f"--endpoint-url={os.getenv('BOTO_ENDPOINT_URL')}")

        command.extend(
            [
                "s3",
                "sync",
                src_dir,
                s3_uri,
            ]
        )
        subprocess.run(
            command,
            check=True,
        )

    def does_object_exist(self, bucket_name: str, object_key: str) -> bool:
        """
        Returns True if the object exists in the bucket and is available to download
        """
        deployment_stage = os.getenv("DEPLOYMENT_STAGE")
        remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX")

        if deployment_stage == "rdev" and remote_dev_prefix is not None:
            # bucket_name is "cellguide-data-public-dev/env-rdev-cellguide"
            # so subpath is "env-rdev-cellguide"
            bucket_name, subpath = bucket_name.split("/", 1)
            object_key = f"{subpath}{remote_dev_prefix}/{object_key}"
        elif deployment_stage == "rdev":
            raise ValueError("REMOTE_DEV_PREFIX must be set when DEPLOYMENT_STAGE is rdev")

        try:
            self.client.head_object(Bucket=bucket_name, Key=object_key)
            return True
        except botocore.exceptions.ClientError as e:
            error_code = int(e.response["Error"]["Code"])
            if error_code == 404:
                return False
            else:
                raise
