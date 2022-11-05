from typing import Optional
from urllib.parse import urlparse

from sqlalchemy.orm import Session

from backend.common.corpora_orm import CollectionVisibility, ProjectLinkType
from backend.common.entities import Collection, Dataset
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    MethodNotAllowedException,
    InvalidParametersHTTPException,
    NotFoundHTTPException,
)
from backend.common.utils.regex import DOI_REGEX_COMPILED
from backend.portal.api.app.v1.authorization import owner_or_allowed


def portal_get_normalized_doi_url(doi_node: dict, errors: list) -> Optional[str]:
    """
    1. Check for DOI uniqueness in the payload
    2. Normalizes it so that the DOI is always a link (starts with https://doi.org)
    3. Returns the newly normalized DOI
    """
    doi_url = doi_node["link_url"]
    parsed = urlparse(doi_url)
    if not parsed.scheme and not parsed.netloc:
        parsed_doi = parsed.path
        if not DOI_REGEX_COMPILED.match(parsed_doi):
            errors.append({"link_type": ProjectLinkType.DOI.name, "reason": "Invalid DOI"})
            return None
        doi_url = f"https://doi.org/{parsed_doi}"
    return doi_url
