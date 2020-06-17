import unittest

from backend.corpora.common.utils.math_utils import sizeof_formatted


class TestMathUtils(unittest.TestCase):
    def test__sizeof_formatted__converts_to_gb(self):
        lots_of_bytes = 4382540267  # Roughly 4.1 GB

        actual_formatted_size = sizeof_formatted(lots_of_bytes)

        expected_formatted_size = "4.1 GiB"
        self.assertEqual(actual_formatted_size, expected_formatted_size)
