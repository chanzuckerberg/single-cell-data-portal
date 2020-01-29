import unittest

from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into


class TestUtilityFunctions(unittest.TestCase):
    def test_merge_two_empty_dictionaries(self):
        dictionary_one = {}
        dictionary_two = {}
        expected_merged_dictionary = {}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_into_an_empty_dictionary(self):
        dictionary_one = {}
        dictionary_two = {"key": "value"}
        expected_merged_dictionary = {"key": ["value"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_empty_dictionary_into_non_empty_dictionary(self):
        dictionary_one = {"key": ["value"]}
        dictionary_two = {}
        expected_merged_dictionary = {"key": ["value"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_two_non_empty_dictionaries_different_keys(self):
        dictionary_one = {"key1": ["value1"]}
        dictionary_two = {"key2": "value2"}
        expected_merged_dictionary = {"key1": ["value1"], "key2": ["value2"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_two_non_empty_dictionaries_same_key(self):
        dictionary_one = {"key1": ["value1"]}
        dictionary_two = {"key1": "value2"}
        expected_merged_dictionary = {"key1": ["value1", "value2"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_two_non_empty_dictionaries_same_and_different_keys(self):
        dictionary_one = {"key1": ["value1"]}
        dictionary_two = {"key1": "value2", "key2": "value3"}
        expected_merged_dictionary = {"key1": ["value1", "value2"], "key2": ["value3"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_two_dictionaries_keys_with_multiple_values(self):
        dictionary_one = {"key1": ["value1"]}
        dictionary_two = {"key1": ["value2", "value3"]}
        expected_merged_dictionary = {"key1": ["value1", "value2", "value3"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)

    def test_merge_empty_dictionary_and_dictionary_with_lists(self):
        dictionary_one = {}
        dictionary_two = {"key1": ["value1", "value2"]}
        expected_merged_dictionary = {"key1": ["value1", "value2"]}

        actual_merged_dictionary = merge_dictionary_into(dictionary_one, dictionary_two)

        self.assertEqual(actual_merged_dictionary, expected_merged_dictionary)
