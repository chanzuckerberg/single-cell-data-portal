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
        data = convert_dataclass_to_dict(data, remove_nones=True)

    with open(path, "w") as f:
        json.dump(data, f)
