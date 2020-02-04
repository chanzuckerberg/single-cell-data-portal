from typing import Dict


def merge_dictionary_into(dict1: Dict, dict2: Dict):
    merged_dictionary = dict1
    for key, value in dict2.items():
        if key in merged_dictionary:
            if type(value) is list:
                merged_dictionary[key] = merged_dictionary[key] + value
            else:
                merged_dictionary[key].append(value)
        else:
            merged_dictionary[key] = value if type(value) is list else [value]
    return merged_dictionary
