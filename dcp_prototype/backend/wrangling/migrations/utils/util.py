import json
import typing
from typing import Dict

import boto3

s3 = boto3.client("s3")


def merge_dictionary_into(dict1: Dict, dict2: Dict):
    merged_dictionary = dict1
    for key, value in dict2.items():
        if key in merged_dictionary:
            if type(value) is list:
                merged_dictionary[key] = merged_dictionary[key] + value
            else:
                merged_dictionary[key].append(value)
        else:
            merged_dictionary[key] = value if type(value) is list else [value]
    return merged_dictionary


def list_files_in_bucket(bucket: str, key: str, continuation_token: str = None, files: typing.Dict[str, str] = {}, checksums: bool=False) -> Dict[str, str]:
    """
    Lists all files in a bucket with the prefix key
    :param bucket: Name of AWS 3s bucket
    :param key: File prefix to check for
    :param continuation_token: s3 will return a max of 1000 files, pass in the continuation token to get the next 1000 for the {bucket}/{prefix}
    :param files: dict of all of the file names or checksums
    :param checksums: Boolean, if true returns list of file checksums
    :return: the file_list, a list of all of the file names or checksums (Depending on value of checksums param)
    """
    if continuation_token:
        response = s3.list_objects_v2(
            Bucket=bucket,
            Prefix=key,
            ContinuationToken=continuation_token
        )
    else:
        response = s3.list_objects_v2(
            Bucket=bucket,
            Prefix=key
        )
    if response['KeyCount'] > 0:
        for x in response['Contents']:
            if checksums:
                checksum = json.loads(x['ETag'])
                files[checksum] = x
            else:

                files[x['Key']] = x
    if response['IsTruncated'] is True:
        list_files_in_bucket(bucket, key, response['NextContinuationToken'], files, checksums)
    return files
