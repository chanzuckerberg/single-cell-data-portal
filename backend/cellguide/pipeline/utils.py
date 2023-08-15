import json
from dataclasses import asdict, is_dataclass


def output_json(data, path):
    """
    data - a dict or a dataclass
    path - str
    """

    if is_dataclass(data):
        data = asdict(data)

    with open(path, "w") as f:
        json.dump(data, f)
