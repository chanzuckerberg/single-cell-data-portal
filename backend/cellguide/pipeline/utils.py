import json
from dataclasses import asdict, is_dataclass


def convert_dataclass_to_dict(data):
    """
    Convert a dataclass or dict of dataclasses to a dict and remove all keys
    with None values.
    """
    if is_dataclass(data):
        data = asdict(data)

    elif isinstance(data, dict):
        for key, value in data.items():
            if is_dataclass(value):
                data[key] = asdict(value)

    return _remove_none_values(data)


def output_json(data, path):
    """
    Output JSON while handling the case where `data` is an instance of a dataclass

    Arguments
    ---------
    data - a dict or a dataclass
    path - str
    """

    if is_dataclass(data):
        data = convert_dataclass_to_dict(data)

    with open(path, "w") as f:
        json.dump(data, f)


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
