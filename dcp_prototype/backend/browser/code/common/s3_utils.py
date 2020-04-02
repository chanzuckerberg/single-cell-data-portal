import os

import boto3
from botocore.exceptions import ClientError


def generate_file_url(file_prefix, expiration=3600, s3=None):
    """
    Generate a presigned URL for a file for user download.
    :param file_prefix: S3 prefix location of the file
    :param expiration: Presigned URL expiration in seconds
    :return: Presigned URL to download the requested file
    """
    if not s3:
        s3 = boto3.client("s3")

    bucket_name = f"dcp-browser-bucket-{os.environ['DEPLOYMENT_STAGE']}"

    try:
        response = s3.generate_presigned_url(
            "get_object", Params={"Bucket": bucket_name, "Key": file_prefix}, ExpiresIn=expiration
        )
    except ClientError as e:
        print(f"Failed to generate presigned URL for {file_prefix} with error: {e}")
        return ""

    return response
