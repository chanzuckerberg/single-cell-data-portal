import json


def compare_dicts(dict1, dict2):
    """
    Recursively compares two dictionaries, handling unordered lists and lists of dictionaries.
    """
    if len(dict1) != len(dict2):
        return False

    for key in dict1:
        if key not in dict2:
            return False

        value1 = dict1[key]
        value2 = dict2[key]

        if isinstance(value1, dict) and isinstance(value2, dict):
            if not compare_dicts(value1, value2):
                return False
        elif isinstance(value1, list) and isinstance(value2, list):
            if len(value1) != len(value2):
                return False

            # Handle lists of dictionaries (sort by keys before comparing)
            if all(isinstance(v, dict) for v in value1) and all(isinstance(v, dict) for v in value2):
                sorted_value1 = sorted(value1, key=lambda d: json.dumps(d, sort_keys=True))
                sorted_value2 = sorted(value2, key=lambda d: json.dumps(d, sort_keys=True))
                if not all(compare_dicts(d1, d2) for d1, d2 in zip(sorted_value1, sorted_value2, strict=False)):
                    return False

            # Handle lists with mixed data (sort elements before comparing)
            elif sorted(value1, key=str) != sorted(value2, key=str):
                return False
        else:
            if value1 != value2:
                return False

    return True


def sort_dataframe(df):
    return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)
