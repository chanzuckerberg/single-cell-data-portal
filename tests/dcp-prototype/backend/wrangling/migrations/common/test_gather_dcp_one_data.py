import random
import unittest

from dcp_prototype.backend.wrangling.migrations.common.gather_dcp_one_data import gather_group_file_list


class TestGatherDcpOneData(unittest.TestCase):
    def test_gather_group_file_list_ensure_json_file_are_first(self):
        json_ending_files = ["a.json", "b.json", "c.json"]
        non_json_ending_file = ["randofile", "someotherfile.tar.gz", "12345-12345-12345-12345"]
        json_and_non_json_files = json_ending_files + non_json_ending_file

        mixed_list_of_files = random.sample(json_and_non_json_files, k=len(json_and_non_json_files))
        ordered_list_of_files = gather_group_file_list(mixed_list_of_files)

        self.assertEqual(set(ordered_list_of_files[0]), set(json_ending_files))
        self.assertEqual(set(ordered_list_of_files[1]), set(non_json_ending_file))

    def test_gather_group_file_list_only_json_files(self):
        json_ending_files = ["a.json", "b.json", "c.json"]

        mixed_list_of_files = random.sample(json_ending_files, k=len(json_ending_files))
        ordered_list_of_files = gather_group_file_list(mixed_list_of_files)

        self.assertEqual(set(ordered_list_of_files[0]), set(json_ending_files))
        self.assertEqual(len(ordered_list_of_files[1]), 0)

    def test_gather_group_file_list_only_non_json_files(self):
        non_json_ending_file = ["randofile", "someotherfile.tar.gz", "12345-12345-12345-12345"]

        mixed_list_of_files = random.sample(non_json_ending_file, k=len(non_json_ending_file))
        ordered_list_of_files = gather_group_file_list(mixed_list_of_files)

        self.assertEqual(len(ordered_list_of_files[0]), 0)
        self.assertEqual(set(ordered_list_of_files[1]), set(non_json_ending_file))
