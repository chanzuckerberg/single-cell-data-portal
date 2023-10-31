import itertools

import pytest

from backend.layers.business.exceptions import InvalidMetadataException
from backend.layers.common.entities import CollectionLinkType, CollectionMetadata, Link
from backend.layers.common.validation import verify_collection_metadata


def test_blank_fields():
    errors = []
    body = CollectionMetadata(name="", contact_name="", description="", contact_email="", links=[])
    with pytest.raises(InvalidMetadataException):
        verify_collection_metadata(body, errors)
    error_message = "Cannot be blank."
    assert {"name": "description", "reason": error_message} in errors
    assert {"name": "name", "reason": error_message} in errors
    assert {"name": "contact_name", "reason": error_message} in errors
    assert {"name": "contact_email", "reason": error_message} in errors


def test_OK():
    body = CollectionMetadata(
        name="something", contact_name="a name", description="description", contact_email="email@place.com", links=[]
    )
    errors = []
    verify_collection_metadata(body, errors)
    assert not errors


@pytest.mark.parametrize(
    "link_type,test_url", itertools.product(CollectionLinkType, ["://", "google", ".com", "google.com", "https://"])
)
def test__link__INVALID(link_type, test_url):
    if link_type.name == "DOI":
        return
    body = CollectionMetadata(
        name="string",
        contact_name="string",
        description="string",
        contact_email="email@mail.com",
        links=[Link(type=link_type.name, uri=test_url, name="")],
    )
    errors = []
    verify_collection_metadata(body, errors)
    expected_error = [dict(reason="Invalid URI.", name="links[0]", value=test_url)]
    assert expected_error == errors


@pytest.mark.parametrize(
    "link_type,test_url",
    itertools.product(CollectionLinkType, ["https://www.google.com", "http://somewhere.org/path/?abcd=123"]),
)
def test__link__OK(link_type, test_url):
    if link_type.name == "DOI":
        return
    body = CollectionMetadata(
        name="string",
        contact_name="string",
        description="string",
        contact_email="email@mail.com",
        links=[Link(type=link_type.name, uri=test_url, name="")],
    )
    errors = []
    verify_collection_metadata(body, errors)
    assert not errors


@pytest.mark.parametrize("bad_email", ["@.", "email@.", "@place.com", "email@.com", "email@place."])
def test_invalid_email(bad_email):
    body = CollectionMetadata(
        name="string", contact_name="string", description="string", contact_email=bad_email, links=[]
    )
    errors = []
    with pytest.raises(InvalidMetadataException):
        verify_collection_metadata(body, errors)
    assert [{"name": "contact_email", "reason": "Invalid format."}] == errors


@pytest.mark.parametrize("invalid_string", [b"\x00some data", b"text\x1f", b"text\x01", b"\x7ftext"])
def test_invalid_characters_in_field(invalid_string):
    errors = []
    string = invalid_string.decode(encoding="utf-8")
    body = CollectionMetadata(
        name=string, contact_name=string, description=string, contact_email="email@email.com", links=[]
    )
    with pytest.raises(InvalidMetadataException):
        verify_collection_metadata(body, errors)
    error_message = "Invalid characters detected."
    assert len(errors) == 1
    assert {"name": "name", "reason": error_message} in errors
