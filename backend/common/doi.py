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


def doi_curie_from_link(doi: str) -> str:
    # Remove the https://doi.org/ (or other) domain part
    parsed = urlparse(doi)
    if parsed.scheme and parsed.netloc:
        doi = parsed.path.lstrip("/")
    return doi


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


def clean_doi(doi: str) -> str:
    """
    Cleans the DOI string. Formats handled:
    - DOI 10.1182/ bloodadvances.2017015073
    - DOI:10.1167/iovs.15-18117
    - DOI: 10.1002/biot.201200199
    - DOI: 10.1111/j.1440-1827.1995.tb03518.x.
    - https://doi.org/10.1101/2021.01.02.425073

    Parameters
    ----------
    doi : str
        The DOI string to be cleaned.

    Returns
    -------
    str
        The cleaned DOI string.
    """
    doi = doi.strip()
    if doi == "No DOI":
        return ""

    # Remove trailing periods from the DOI. This handles the
    # "10.1111/j.1440-1827.1995.tb03518.x."-type cases.
    if doi != "" and doi[-1] == ".":
        doi = doi[:-1]

    # Remove any invalid tokens from the DOI. Invalid tokens include:
    # "DOI", "DOI:", "DOI: ", and "https://doi.org/".
    regex = re.compile(r"\bDOI[: ]?\s*|https://doi.org/", re.IGNORECASE)
    doi = regex.sub("", doi)

    # Remove all whitespace from the DOI. This handles the
    # "10.1182/ bloodadvances.2017015073"-type cases, as well as any other
    # leading or trailing whitespace.
    doi = re.sub(r"\s+", "", doi.strip())

    return doi
