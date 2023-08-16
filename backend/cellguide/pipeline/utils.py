import json
from dataclasses import is_dataclass

from backend.common.utils.dataclass import convert_dataclass_to_dict


def output_json(data, path):
    """
    Output JSON while handling the case where `data` is an instance of a dataclass

    Arguments
    ---------
    data - a dict or a dataclass
    path - str
    """

    if is_dataclass(data):
        data = convert_dataclass_to_dict_and_strip_nones(data)

    with open(path, "w") as f:
        json.dump(data, f)


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
