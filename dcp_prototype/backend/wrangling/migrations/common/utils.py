import os
from urllib.parse import urlparse


def set_attribute_value(entity, attr, value):
    """
    Given an entity object, set its member variable "attr" to the given value.
    """

    setattr(entity, attr, str(value))


def append_value_to_attribute(entity, attr, value):
    """
    Given an entity object instance, append "value" to its member variable "attr." Attr is expected to be a list of
    values of the same type as the given value object.
    """

    if type(getattr(entity, attr)) is not list:
        raise RuntimeError(
            f"ERROR: Attempted to append a value {value} to a non-list attribute {attr} for entity object {entity}."
        )

    getattr(entity, attr).append(value)


def append_unique_value_to_attribute(entity, attr, value):
    """
    Given an entity object instance, append "value" to its member variable "attr" only if it does not already exist.
    Attr is expected to be a list of values of the same type as the given value object.
    """

    if type(getattr(entity, attr)) is not list:
        raise RuntimeError(
            f"ERROR: Attempted to append a value {value} to a non-list attribute {attr} for entity object {entity}."
        )

    if value not in getattr(entity, attr):
        getattr(entity, attr).append(value)


def get_entity_type(object_name, object_data):
    """
    Given a dictionary representation of an instance of a metadata entity, extract and return its entity type.
    """

    described_by = object_data.get("describedBy")
    if described_by:
        path = urlparse(described_by).path
        entity_type = os.path.basename(path)
        return entity_type
    else:
        raise RuntimeError(f"ERROR: Object {object_name} does not have a describedBy field.")
