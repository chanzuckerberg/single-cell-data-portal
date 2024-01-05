import re
from typing import Optional
from urllib.parse import urlparse

from backend.common.utils.regex import CURIE_REFERENCE_REGEX, DOI_REGEX_COMPILED


def get_doi_link_node(body: dict, errors: list) -> Optional[dict]:
    links = body.get("links", [])
    dois = [link for link in links if link["link_type"] == "DOI"]

    if not dois:
        return None

    # Verify that a single DOI exists
    if len(dois) > 1:
        errors.append({"link_type": "DOI", "reason": "Can only specify a single DOI"})
        return None

    return dois[0]


def curation_get_normalized_doi_url(doi: str, errors: list) -> Optional[str]:
    # Regex below is adapted from https://bioregistry.io/registry/doi 'Pattern for CURIES'
    if not re.match(CURIE_REFERENCE_REGEX, doi):
        errors.append({"name": "DOI", "reason": "DOI must be a CURIE reference."})
        return None
    return f"https://doi.org/{doi}"


def portal_get_normalized_doi_url(doi_node: dict, errors: list) -> Optional[str]:
    """
    1. Check for DOI uniqueness in the payload
    2. Normalizes it so that the DOI is always a link (starts with https://doi.org)
    3. Returns the newly normalized DOI
    """
    doi_url = doi_node["link_url"]
    parsed = urlparse(doi_url)
    if not parsed.scheme and not parsed.netloc:
        parsed_doi = parsed.path.lstrip("/")
        if not DOI_REGEX_COMPILED.match(parsed_doi):
            errors.append({"link_type": "DOI", "reason": "Invalid DOI"})
            return None
        doi_url = f"https://doi.org/{parsed_doi}"
    return doi_url
