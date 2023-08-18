import copy
from dataclasses import asdict, is_dataclass


def convert_dataclass_to_dict(data):
    """
    Convert a dataclass or dict of dataclasses to a dict.

    Arguments
    ---------
    data - dataclass or dict of dataclasses
    """
    data_copy = copy.deepcopy(data)

    if is_dataclass(data_copy):
        data_copy = asdict(data_copy)

    elif isinstance(data_copy, dict):
        for key, value in data_copy.items():
            if is_dataclass(value):
                data_copy[key] = asdict(value)
            elif isinstance(value, list):
                data_copy[key] = [asdict(i) if is_dataclass(i) else i for i in value]

    return data_copy
