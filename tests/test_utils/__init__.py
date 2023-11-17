import contextlib
from contextlib import contextmanager


class CompareDictsAddin:
    def assert_dicts_equal(self, dict1, dict2, fail_fast=True, max_diff=None):
        """
        This function recursive compares two dictionaries handling cases where the values
        are unordered arrays with elements that could be dictionaries.
        :param dict1: first dictionary
        :param dict2: second dictionary
        :param fail_fast: if True, raise an exception on the first mismatch
        :param max_diff: maximum number of differences to report
        """
        _mismatch = []
        _path = ["$"]

        class Stop(Exception):
            pass

        def update_mismatch(message):
            _mismatch.append(message)
            if fail_fast:
                raise Stop()

        @contextmanager
        def path(key):
            try:
                yield _path.append(key)
            finally:
                _path.pop()

        def _compare_dict(_dict1, _dict2):
            if len(_dict1) != len(_dict2):
                update_mismatch("{}: Lengths differ: {} vs {}".format(".".join(path), len(_dict1), len(_dict2)))

            for key in dict1:
                if key not in _dict2:
                    update_mismatch("{}: Key {} not in dict2".format(".".join(path), key))
                    continue
                value1 = _dict1[key]
                value2 = _dict2[key]

                with path(key):
                    if isinstance(value1, dict) and isinstance(value2, dict):
                        _compare_dict(value1, value2)
                    elif isinstance(value1, list) and isinstance(value2, list):
                        if len(value1) != len(value2):
                            update_mismatch(
                                "{}: Lengths differ: {} vs {}".format(".".join(path), len(value1), len(value2))
                            )
                        # check if the lists contain dictionaries as elements
                        elif len(value1) > 0 and isinstance(value1[0], dict) and isinstance(value2[0], dict):
                            for i in range(len(value1)):
                                with path(f"[{i}]"):
                                    _compare_dict(value1[i], value2[i])
                        elif sorted(value1) != sorted(value2):
                            # compare as sets and print the differences
                            value1_set = set(value1)
                            value2_set = set(value2)
                            update_mismatch(
                                "{}: Values differ: {} vs {}".format(
                                    ".".join(path), value1_set - value2_set, value2_set - value1_set
                                )
                            )
                    else:
                        if value1 != value2:
                            update_mismatch("{}: Values differ: {} vs {}".format(".".join(path), value1, value2))

        with contextlib.suppress(Stop):
            _compare_dict(dict1, dict2)

        if _mismatch:
            if len(_mismatch) > max_diff:
                raise self.failureException(
                    f"Mismatched keys: {len[_mismatch]}" + "\n\t".join(_mismatch[:max_diff]) + "\n\t..."
                )
            else:
                raise self.failureException(f"Mismatched keys: {len[_mismatch]}" + "\n\t".join(_mismatch))


def sort_dataframe(df):
    return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)
