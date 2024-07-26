import re
from typing import List, Union
from urllib.parse import urlparse

from backend.layers.business.entities import CollectionMetadataUpdate
from backend.layers.business.exceptions import InvalidMetadataException
from backend.layers.common.entities import CollectionMetadata, DatasetArtifactMetadataUpdate, Link
from backend.layers.common.regex import CONTROL_CHARS, EMAIL_REGEX

control_char_re = re.compile(CONTROL_CHARS)

valid_consortia = {
    "Allen Institute for Brain Science",
    "BRAIN Initiative",
    "CZ Biohub",
    "CZI Neurodegeneration Challenge Network",
    "CZI Single-Cell Biology",
    "European Union’s Horizon 2020",
    "GenitoUrinary Development Molecular Anatomy Project (GUDMAP)",
    "Gut Cell Atlas",
    "Human BioMolecular Atlas Program (HuBMAP)",
    "Human Cell Atlas (HCA)",
    "Human Pancreas Analysis Program (HPAP)",
    "Human Tumor Atlas Network (HTAN)",
    "Kidney Precision Medicine Project (KPMP)",
    "LungMAP",
    "Pediatric Center of Excellence in Nephrology (PCEN)",
    "SEA-AD",
    "Wellcome HCA Strategic Science Support",
}


def _verify_field(
    metadata: Union[CollectionMetadata, CollectionMetadataUpdate, DatasetArtifactMetadataUpdate],
    key: str,
    check_existence: bool,
    errors: list,
) -> None:
    """
    Verifies field exists and is not empty.

    :param metadata: Metadata update to validate.
    :param key: Field to validate in metadata update.
    :param check_existence: If True, field must exist in metadata update.
    :param errors: List of errors to append to.
    """
    value = getattr(metadata, key)
    if check_existence and value is None:
        # if checks_existence is true, value cannot be None since it must be required
        errors.append({"name": key, "reason": "Cannot be empty."})
    elif value is not None and not value:
        # In any case, if a value is defined, it cannot be falsey (aka blank)
        errors.append({"name": key, "reason": "Cannot be blank."})
    elif value is not None and (key == "name" or key == "title") and control_char_re.search(value):
        errors.append({"name": key, "reason": "Invalid characters detected."})
    else:
        return value


def _verify_collection_metadata_fields(
    metadata: Union[CollectionMetadata, CollectionMetadataUpdate], check_existence: bool, errors: list
) -> None:
    """
    Verifies if the fields in CollectionMetadata or CollectionMetadataUpdate are valid. They must:
    - Be non empty
    - Not contain invalid characters
    - If the field is an email, it should be in the right format
    """

    def verify_collection_consortia(metadata: Union[CollectionMetadata, CollectionMetadataUpdate], errors: list):
        consortia = metadata.consortia
        if consortia:
            for consortium in consortia:
                if consortium not in valid_consortia:
                    errors.append({"name": "consortia", "reason": "Invalid consortia."})

    contact_email = _verify_field(metadata, "contact_email", check_existence, errors)
    if contact_email:
        result = EMAIL_REGEX.match(contact_email)
        if not result:
            errors.append({"name": "contact_email", "reason": "Invalid format."})

    _verify_field(metadata, "description", check_existence, errors)
    _verify_field(metadata, "name", check_existence, errors)
    _verify_field(metadata, "contact_name", check_existence, errors)

    verify_collection_consortia(metadata, errors)


def verify_collection_links(links: List[Link], errors: list) -> None:
    def _error_message(i: int, _url: str) -> dict:
        return {"name": f"links[{i}]", "reason": "Invalid URL.", "value": _url}

    for index, link in enumerate(links):
        url = link.uri
        try:
            result = urlparse(url.strip())
            if not all([result.scheme, result.netloc]):
                errors.append(_error_message(index, url))
        except ValueError:
            errors.append(_error_message(index, url))


def verify_collection_metadata_update(metadata: CollectionMetadataUpdate, errors: list) -> None:
    _verify_collection_metadata_fields(metadata, check_existence=False, errors=errors)
    if errors:
        raise InvalidMetadataException(errors=errors)
    if metadata.links is not None:
        verify_collection_links(metadata.links, errors)


def verify_collection_metadata(metadata: CollectionMetadata, errors: list) -> None:
    """
    Verify if `CollectionMetadata` is well defined. Since `CollectionMetadataCreate
    """
    _verify_collection_metadata_fields(metadata, check_existence=True, errors=errors)
    if errors:
        raise InvalidMetadataException(errors=errors)
    verify_collection_links(metadata.links, errors)


def verify_dataset_artifact_metadata_update(metadata_update: DatasetArtifactMetadataUpdate) -> None:
    """
    Verify values of `DatasetArtifactMetadataUpdate` are valid. Currently only title is available for
    update via the FE and DP and Discover APIs; title must be specified, and must not contain special characters.

    :param metadata_update: Metadata update to validate.
    """
    errors = []
    _verify_field(metadata_update, "title", True, errors)
    if errors:
        raise InvalidMetadataException(errors=errors)
