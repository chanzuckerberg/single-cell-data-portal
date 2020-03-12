import unittest

from dcp_prototype.backend.wrangling.migrations.common.utils import set_attribute_value, get_entity_type, \
    append_value_to_attribute, append_unique_value_to_attribute


class TestUtils(unittest.TestCase):
    def test_set_attribute_value_successful(self):
        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = None

        dummy_instance = DummyClass()
        dummy_variable = "dummy_variable"
        dummy_value = "dummy_value"

        set_attribute_value(dummy_instance, dummy_variable, dummy_value)

        self.assertEqual(dummy_instance.dummy_variable, dummy_value)

    def test_set_attribute_value_non_existent_variable_successful(self):
        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = None

        dummy_instance = DummyClass()
        dummy_variable = "non_existent_dummy_variable"
        dummy_value = "dummy_value"

        set_attribute_value(dummy_instance, dummy_variable, dummy_value)

        self.assertEqual(dummy_instance.non_existent_dummy_variable, dummy_value)

    def test_append_value_to_attribute_successful(self):
        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = []

        dummy_instance = DummyClass()
        dummy_variable = "dummy_variable"
        dummy_value = "a cool value"

        append_value_to_attribute(dummy_instance, dummy_variable, dummy_value)

        expected_value = [dummy_value]
        self.assertEqual(dummy_instance.dummy_variable, expected_value)

    def test_append_value_to_attribute_with_existing_values_successful(self):
        dummy_variable = "dummy_variable"
        dummy_value_existing = "some cool value"
        dummy_value_new = "another cool value"

        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = [dummy_value_existing]

        dummy_instance = DummyClass()
        append_value_to_attribute(dummy_instance, dummy_variable, dummy_value_new)

        expected_value = [dummy_value_existing, dummy_value_new]
        self.assertEqual(dummy_instance.dummy_variable, expected_value)

    def test_append_value_to_attribute_non_list_attribute_fails(self):
        dummy_variable = "dummy_variable"
        dummy_value_existing = "some cool value"
        dummy_value_new = "another cool value"

        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = dummy_value_existing

        dummy_instance = DummyClass()

        with self.assertRaises(RuntimeError) as error:
            append_value_to_attribute(dummy_instance, dummy_variable, dummy_value_new)

            self.assertIn(f"Attempted to append a value {dummy_value_new} to a non-list attribute", error)

    def test_append_unique_value_to_attribute_successful(self):
        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = []

        dummy_instance = DummyClass()
        dummy_variable = "dummy_variable"
        dummy_value = "a cool value"

        append_unique_value_to_attribute(dummy_instance, dummy_variable, dummy_value)

        expected_value = [dummy_value]
        self.assertEqual(dummy_instance.dummy_variable, expected_value)

    def test_append_unique_value_to_attribute_with_preexisting_successful(self):
        dummy_variable = "dummy_variable"
        dummy_value_existing = "some cool value"
        dummy_value_new = "another cool value"

        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = [dummy_value_existing]

        dummy_instance = DummyClass()
        append_unique_value_to_attribute(dummy_instance, dummy_variable, dummy_value_new)

        expected_value = [dummy_value_existing, dummy_value_new]
        self.assertEqual(dummy_instance.dummy_variable, expected_value)

    def test_append_unique_value_to_attribute_with_preexisting_same_value_successful(self):
        dummy_variable = "dummy_variable"
        dummy_value = "another cool value"

        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = [dummy_value]

        dummy_instance = DummyClass()
        append_unique_value_to_attribute(dummy_instance, dummy_variable, dummy_value)

        expected_value = [dummy_value]
        self.assertEqual(dummy_instance.dummy_variable, expected_value)

    def test_append_unique_value_to_attribute_non_list_attribute_fails(self):
        dummy_variable = "dummy_variable"
        dummy_value_existing = "some cool value"
        dummy_value_new = "another cool value"

        # Create a dummy class
        class DummyClass:
            def __init__(self):
                self.dummy_variable = dummy_value_existing

        dummy_instance = DummyClass()

        with self.assertRaises(RuntimeError) as error:
            append_unique_value_to_attribute(dummy_instance, dummy_variable, dummy_value_new)

            self.assertIn(f"Attempted to append a value {dummy_value_new} to a non-list attribute", error)

    def test_get_entity_type_successful(self):
        expected_entity_type = "donor_organism"
        entity_name = "donor_organism_12345.json"
        entity_data = {
            "describedBy": f"https://schema.humancellatlas.org/type/biomaterial/15.5.0/{expected_entity_type}",
            "another_exciting_field": "some exciting value"
        }

        actual_entity_type = get_entity_type(entity_name, entity_data)

        self.assertEqual(actual_entity_type, expected_entity_type)

    def test_get_entity_type_no_described_by_field_fails(self):
        entity_name = "donor_organism_12345.json"
        entity_data = {
            "another_exciting_field": "some exciting value"
        }

        with self.assertRaises(RuntimeError) as error:
            get_entity_type(entity_name, entity_data)

            self.assertIn("does not have a describedBy field", error)
