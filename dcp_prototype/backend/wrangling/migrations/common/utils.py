import json
import os
import typing
import boto3

from urllib.parse import urlparse


def set_attribute_value(entity, attr, value):
    """
    Given an entity object, set its member variable "attr" to the given value.
    """

    setattr(entity, attr, str(value))


def append_value_to_attribute(entity, attr, value):
    """
    Given an entity object instance, append "value" to its member variable "attr." Attr is expected to be a list of
    values of the same type as the given value object.
    """

    if type(getattr(entity, attr)) is not list:
        raise RuntimeError(
            f"ERROR: Attempted to append a value {value} to a non-list attribute {attr} for entity object {entity}."
        )

    getattr(entity, attr).append(value)


def append_unique_value_to_attribute(entity, attr, value):
    """
    Given an entity object instance, append "value" to its member variable "attr" only if it does not already exist.
    Attr is expected to be a list of values of the same type as the given value object.
    """

    if type(getattr(entity, attr)) is not list:
        raise RuntimeError(
            f"ERROR: Attempted to append a value {value} to a non-list attribute {attr} for entity object {entity}."
        )

    if value not in getattr(entity, attr):
        getattr(entity, attr).append(value)


def get_entity_type(object_name, object_data):
    """
    Given a dictionary representation of an instance of a metadata entity, extract and return its entity type.
    """

    described_by = object_data.get("describedBy")
    if described_by:
        path = urlparse(described_by).path
        entity_type = os.path.basename(path)
        return entity_type
    else:
        raise RuntimeError(f"ERROR: Object {object_name} does not have a describedBy field.")


def merge_dictionary_into(dict1: typing.Dict, dict2: typing.Dict):
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


def list_files_in_bucket(bucket, key, continuation_token=None, files={}, checksums=False):
    """
    Lists all files in a bucket with the prefix key
    :param bucket: Name of AWS S3 bucket
    :param key: File prefix to check for
    :param continuation_token: s3 will return a max of 1000 files, pass in the continuation token to get the next 1000
    for the {bucket}/{prefix}
    :param files: dict of all of the file names or checksums
    :param checksums: Boolean, if true returns list of file checksums
    :return: the file_list, a list of all of the file names or checksums (Depending on value of checksums param)
    """
    s3 = boto3.client("s3")

    if continuation_token:
        response = s3.list_objects_v2(Bucket=bucket, Prefix=key, ContinuationToken=continuation_token)
    else:
        response = s3.list_objects_v2(Bucket=bucket, Prefix=key)
    if response["KeyCount"] > 0:
        for x in response["Contents"]:
            if checksums:
                checksum = json.loads(x["ETag"])
                files[checksum] = x
            else:

                files[x["Key"]] = x
    if response["IsTruncated"] is True:
        # due to python recursion limit (1000) this may fail when more than 1 million data files have the
        # same prefix (key)
        list_files_in_bucket(bucket, key, response["NextContinuationToken"], files, checksums)
    return files
