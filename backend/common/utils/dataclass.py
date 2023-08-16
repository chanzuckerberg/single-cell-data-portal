from dataclasses import asdict, is_dataclass


def convert_dataclass_to_dict(data):
    """
    Convert a dataclass or dict of dataclasses to a dict.

    Arguments
    ---------
    data - dataclass or dict of dataclasses
    """
    if is_dataclass(data):
        data = asdict(data)

    elif isinstance(data, dict):
        for key, value in data.items():
            if is_dataclass(value):
                data[key] = asdict(value)

    return data
