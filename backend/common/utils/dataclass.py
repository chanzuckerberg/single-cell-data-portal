from dataclasses import asdict, is_dataclass


def convert_dataclass_to_dict(data, remove_nones=True):
    """
    Convert a dataclass or dict of dataclasses to a dict and optionally
    remove all keys with None values (default behavior).

    Arguments
    ---------
    data - dataclass or dict of dataclasses
    remove_nones - bool, if True, remove all keys with None values
    """
    if is_dataclass(data):
        data = asdict(data)

    elif isinstance(data, dict):
        for key, value in data.items():
            if is_dataclass(value):
                data[key] = asdict(value)

    return _remove_none_values(data)


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
