import re

from urllib.parse import urlparse

from backend.layers.common.entities import CollectionLinkType, CollectionMetadata
from backend.layers.common.regex import CONTROL_CHARS, EMAIL_REGEX

control_char_re = re.compile(CONTROL_CHARS)

def verify_collection_metadata(metadata: CollectionMetadata, errors: list) -> None:
    def check(key):
        if key in body.keys():
            if not body[key]:
                errors.append({"name": key, "reason": "Cannot be blank."})
            elif key == "name" and control_char_re.search(body[key]):
                errors.append({"name": key, "reason": "Invalid characters detected."})
            else:
                return body[key]

    contact_email = check("contact_email")
    if contact_email:
        result = EMAIL_REGEX.match(contact_email)
        if not result:
            errors.append({"name": "contact_email", "reason": "Invalid format."})

    check("description")
    check("name")
    check("contact_name")

    verify_collection_links(body, errors)

def verify_collection_links(body: dict, errors: list) -> None:
    def _error_message(i: int, _url: str) -> dict:
        return {"name": f"links[{i}]", "reason": "Invalid URL.", "value": _url}

    for index, link in enumerate(body.get("links", [])):
        if link["link_type"] == CollectionLinkType.DOI:
            continue
        url = link["link_url"]
        try:
            result = urlparse(url.strip())
        except ValueError:
            errors.append(_error_message(index, url))
        if not all([result.scheme, result.netloc]):
            errors.append(_error_message(index, url))