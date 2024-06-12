import gzip
import json
import logging
import os

from backend.common.utils.dataclass import convert_dataclass_to_dict

logger = logging.getLogger(__name__)


def output_json(data, path):
    """
    Output JSON while handling the case where `data` is an instance of a dataclass

    Arguments
    ---------
    data - a dict or a dataclass
    path - str
    """

    data = convert_dataclass_to_dict_and_strip_nones(data)

    os.makedirs(os.path.dirname(path), exist_ok=True)

    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            json.dump(data, f)
    else:
        with open(path, "w") as f:
            json.dump(data, f)


def output_json_per_key(data, path):
    """
    Output each value in the input data to a separate JSON file.
    """

    data = convert_dataclass_to_dict_and_strip_nones(data)

    os.makedirs(os.path.dirname(f"{path}/"), exist_ok=True)

    for key in data:
        filename = key.replace(":", "_")
        output_json(data[key], f"{path}/{filename}.json")


def convert_dataclass_to_dict_and_strip_nones(data):
    """
    A wrapper for dataclass-to-dict conversion that deletes keys with None as values
    """
    return _remove_none_values(convert_dataclass_to_dict(data))


def _remove_none_values(data):
    """
    Recursively remove keys with None values.

    Arguments
    ---------
    data - dict
    """
    if isinstance(data, dict):
        return {k: _remove_none_values(v) for k, v in data.items() if v is not None}
    elif isinstance(data, list):
        return [_remove_none_values(i) for i in data]
    else:
        return data
