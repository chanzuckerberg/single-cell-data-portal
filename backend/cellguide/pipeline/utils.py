import json
from dataclasses import asdict, is_dataclass


def remove_none_values(data):
    if isinstance(data, dict):
        return {k: remove_none_values(v) for k, v in data.items() if v is not None}
    elif isinstance(data, list):
        return [remove_none_values(i) for i in data]
    else:
        return data


def output_json(data, path):
    """
    data - a dict or a dataclass
    path - str
    """

    if is_dataclass(data):
        data = asdict(data)
        # remove values that are set to None as these were unset in the dataclass
        data = remove_none_values(data)

    with open(path, "w") as f:
        json.dump(data, f)
