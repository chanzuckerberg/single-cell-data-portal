import os

import boto3
from botocore.exceptions import ClientError


def generate_file_url(file_prefix: str, expiration: int = 3600, s3=None) -> str:
    """
    Generate a presigned URL for a file for user download.
    :param file_prefix: S3 prefix location of the file
    :param expiration: Presigned URL expiration in seconds
    :param s3: Optional S3 client instance
    :return: Presigned URL to download the requested file
    """
    if not s3:
        s3 = boto3.client("s3")

    bucket_name = f"corpora-data-{os.environ['DEPLOYMENT_STAGE']}"

    try:
        response = s3.generate_presigned_url(
            "get_object", Params={"Bucket": bucket_name, "Key": file_prefix}, ExpiresIn=expiration
        )
    except ClientError as e:
        print(f"Failed to generate presigned URL for {file_prefix} with error: {e}")
        return ""

    return response

def head_file(file_prefix: str, s3=None):
    """
    The HEAD operation retrieves metadata from an object.
    :param file_prefix: S3 prefix location of the file
    :param s3: Optional S3 client instance
    :return: Presigned URL to download the requested file
    """
    if not s3:
        s3 = boto3.client("s3")
    bucket_name = f"corpora-data-{os.environ['DEPLOYMENT_STAGE']}"
    try:
        response = s3.head_object(Bucket=bucket_name, Key=file_prefix)
    except ClientError as e:
        print(f"Failed for {file_prefix} with error: {e}")
        return ""

    return response
