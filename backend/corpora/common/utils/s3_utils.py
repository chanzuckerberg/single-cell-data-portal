import logging
import typing

import boto3
from botocore.exceptions import ClientError

logger = logging.getLogger(__name__)


def generate_file_url(bucket_name: str, file_prefix: str, expiration: int = 3600) -> typing.Union[str, None]:
    """
    Generate a presigned URL for a file for user download.
    :param bucket_name: S3 name of the bucket
    :param file_prefix: S3 prefix location of the file
    :param expiration: Presigned URL expiration in seconds
    :return: Presigned URL to download the requested file
    """
    s3 = boto3.client("s3")

    try:
        response = s3.generate_presigned_url(
            "get_object", Params={"Bucket": bucket_name, "Key": file_prefix}, ExpiresIn=expiration
        )
    except ClientError:
        logger.exception(f"Failed to generate presigned URL for '{file_prefix}'.")
        return

    return response


def head_file(bucket_name: str, file_prefix: str) -> typing.Union[dict, None]:
    """
    The HEAD operation retrieves metadata from an S3 object.
    :param bucket_name: S3 name of the bucket
    :param file_prefix: S3 prefix location of the file
    :return: S3 metadata for the file
    """
    s3 = boto3.client("s3")

    try:
        response = s3.head_object(Bucket=bucket_name, Key=file_prefix)
    except ClientError:
        logger.exception(f"Failed to retrieve meta data for '{file_prefix}'.")
        return

    return response
