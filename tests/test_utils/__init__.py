class CompareDictsAddin:
    _maxDiff = 10

    def assert_dicts_equal(self, dict1, dict2):
        """
        This function recursive compares two dictionaries handling cases where the values
        are unordered arrays with elements that could be dictionaries.
        """
        mismatch = []
        path = ["$"]

        def _compare_dict(_dict1, _dict2):
            if len(_dict1) != len(_dict2):
                mismatch.append("{}: Lengths differ: {} vs {}".format(".".join(path), len(_dict1), len(_dict2)))

            for key in dict1:
                if key not in _dict2:
                    mismatch.append("{}: Key {} not in dict2".format(".".join(path), key))
                    continue
                path.append(key)
                value1 = _dict1[key]
                value2 = _dict2[key]

                if isinstance(value1, dict) and isinstance(value2, dict):
                    _compare_dict(value1, value2)
                elif isinstance(value1, list) and isinstance(value2, list):
                    if len(value1) != len(value2):
                        mismatch.append("{}: Lengths differ: {} vs {}".format(".".join(path), len(value1), len(value2)))
                    # check if the lists contain dictionaries as elements
                    elif len(value1) > 0 and isinstance(value1[0], dict) and isinstance(value2[0], dict):
                        for i in range(len(value1)):
                            path.append(f"[{i}]")
                            _compare_dict(value1[i], value2[i])
                            path.pop()
                    elif sorted(value1) != sorted(value2):
                        # compare as sets and print the differences
                        value1_set = set(value1)
                        value2_set = set(value2)
                        mismatch.append(
                            "{}: Values differ: {} vs {}".format(
                                ".".join(path), value1_set - value2_set, value2_set - value1_set
                            )
                        )
                else:
                    if value1 != value2:
                        mismatch.append("{}: Values differ: {} vs {}".format(".".join(path), value1, value2))
                path.pop()

        _compare_dict(dict1, dict2)
        if mismatch:
            if len(mismatch) > self._maxDiff:
                raise self.failureException(
                    f"Mismatched keys: {len[mismatch]}" + "\n\t".join(mismatch[: self._maxDiff]) + "\n\t..."
                )
            else:
                raise self.failureException(f"Mismatched keys: {len[mismatch]}" + "\n\t".join(mismatch))


def sort_dataframe(df):
    return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)
