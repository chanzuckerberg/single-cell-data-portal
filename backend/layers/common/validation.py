import re
from typing import List

from urllib.parse import urlparse

from backend.layers.common.entities import CollectionLinkType, CollectionMetadata, Link
from backend.layers.common.regex import CONTROL_CHARS, EMAIL_REGEX

control_char_re = re.compile(CONTROL_CHARS)


# TODO: this method should also verify CollectionMetadataUpdate
def verify_collection_metadata(metadata: CollectionMetadata, errors: list) -> None:
    def check(key):
        value = getattr(metadata, key)
        if not value:
            errors.append({"name": key, "reason": "Cannot be blank."})
        elif key == "name" and control_char_re.search(value):
            errors.append({"name": key, "reason": "Invalid characters detected."})
        else:
            return value

    contact_email = check("contact_email")
    if contact_email:
        result = EMAIL_REGEX.match(contact_email)
        if not result:
            errors.append({"name": "contact_email", "reason": "Invalid format."})

    check("description")
    check("name")
    check("contact_name")

    verify_collection_links(metadata.links, errors)


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
