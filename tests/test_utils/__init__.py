def compare_dicts(dict1, dict2):
    """
    This function recursive compares two dictionaries handling cases where the values
    are unordered arrays with elements that could be dictionaries.
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
            # check if the lists contain dictionaries as elements
            if len(value1) > 0 and isinstance(value1[0], dict) and isinstance(value2[0], dict):
                for i in range(len(value1)):
                    if not compare_dicts(value1[i], value2[i]):
                        return False
            elif sorted(value1) != sorted(value2):
                return False
        else:
            if value1 != value2:
                return False

    return True


def sort_dataframe(df):
    return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)
