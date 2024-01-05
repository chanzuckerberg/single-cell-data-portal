import html
import logging
from datetime import datetime
from typing import Optional, Tuple
from urllib.parse import urlparse

import requests

from backend.common.corpora_config import CorporaConfig


class CrossrefProviderInterface:
    def doi_curie_from_link(self, doi: str) -> str:
        pass

    def fetch_metadata(self, doi: str) -> Tuple[Optional[dict], Optional[str]]:
        return None, None

    def fetch_preprint_published_doi(self, doi):
        pass


class CrossrefException(Exception):
    pass


class CrossrefFetchException(CrossrefException):
    pass


class CrossrefDOINotFoundException(CrossrefException):
    pass


class CrossrefParseException(CrossrefException):
    pass


class CrossrefProvider(CrossrefProviderInterface):
    """
    Provider class used to call Crossref and retrieve publisher metadata
    """

    def __init__(self) -> None:
        self.base_crossref_uri = "https://api.crossref.org/works"
        try:
            self.crossref_api_key = CorporaConfig().crossref_api_key
        except RuntimeError:
            self.crossref_api_key = None
        super().__init__()

    @staticmethod
    def parse_date_parts(obj):
        date_parts = obj["date-parts"][0]
        year = date_parts[0]
        month = date_parts[1] if len(date_parts) > 1 else 1
        day = date_parts[2] if len(date_parts) > 2 else 1
        return (year, month, day)

    def _fetch_crossref_payload(self, doi):

        if self.crossref_api_key is None:
            logging.info("No Crossref API key found, skipping metadata fetching.")
            return None

        try:
            res = requests.get(
                f"{self.base_crossref_uri}/{doi}",
                headers={"Crossref-Plus-API-Token": f"Bearer {self.crossref_api_key}"},
            )
            res.raise_for_status()
        except requests.RequestException as e:
            if e.response is not None and e.response.status_code == 404:
                raise CrossrefDOINotFoundException from e
            else:
                raise CrossrefFetchException("Cannot fetch metadata from Crossref") from e

        return res

    def doi_curie_from_link(self, doi: str) -> str:
        # Remove the https://doi.org part
        parsed = urlparse(doi)
        if parsed.scheme and parsed.netloc:
            doi = parsed.path.lstrip("/")
        return doi

    def fetch_metadata(self, doi: str) -> Tuple[Optional[dict], Optional[str]]:
        """
        Fetches and extracts publisher metadata from Crossref for a specified DOI.
        If the Crossref API URI isn't in the configuration, we will just return an empty object.
        This is to avoid calling Crossref in non-production environments.
        :param doi: str - DOI uri link or curie identifier
        return: tuple - publisher metadata dict and DOI curie identifier
        """
        doi_curie = self.doi_curie_from_link(doi)

        res = self._fetch_crossref_payload(doi_curie)
        if not res:
            return None, None

        try:
            message = res.json()["message"]

            # Date
            published_date = (
                message.get("published-print") or message.get("published") or message.get("published-online")
            )

            if published_date is None:
                raise CrossrefParseException("Date node missing")

            published_year, published_month, published_day = self.parse_date_parts(published_date)

            dates = []
            for k, v in message.items():
                if isinstance(v, dict) and "date-parts" in v:
                    dt = v["date-parts"][0]
                    dates.append(f"{k}: {dt}")

            # Journal
            try:
                if "short-container-title" in message and message["short-container-title"]:
                    raw_journal = message["short-container-title"][0]
                elif "container-title" in message and message["container-title"]:
                    raw_journal = message["container-title"][0]
                elif "institution" in message:
                    raw_journal = message["institution"][0]["name"]
            except Exception:
                raise CrossrefParseException("Journal node missing") from None

            journal = html.unescape(raw_journal)

            # Authors
            # Note: make sure that the order is preserved, as it is a relevant information
            authors = message["author"]
            parsed_authors = []
            for author in authors:
                if "given" in author and "family" in author:
                    parsed_authors.append({"given": author["given"], "family": author["family"]})
                elif "family" in author:
                    # Assume family is consortium
                    parsed_authors.append({"name": author["family"]})
                elif "name" in author:
                    parsed_authors.append({"name": author["name"]})

            # Preprint
            is_preprint = message.get("subtype") == "preprint"
            if is_preprint:
                published_metadata, published_doi_curie = self.fetch_published_metadata(message)
                if published_metadata and published_doi_curie:  # if not, use preprint doi curie
                    return published_metadata, published_doi_curie

            return (
                {
                    "authors": parsed_authors,
                    "published_year": published_year,
                    "published_month": published_month,
                    "published_day": published_day,
                    "published_at": datetime.timestamp(datetime(published_year, published_month, published_day)),
                    "journal": journal,
                    "is_preprint": is_preprint,
                },
                doi_curie,
            )
        except Exception as e:
            raise CrossrefParseException("Cannot parse metadata from Crossref") from e

    def fetch_published_metadata(self, doi_response_message: dict) -> Tuple[Optional[dict], Optional[str]]:
        try:
            published_doi = doi_response_message["relation"]["is-preprint-of"]
            # the new DOI to query for ...
            for entity in published_doi:
                if entity["id-type"] == "doi":
                    return self.fetch_metadata(entity["id"])
        except Exception:  # if fetch of published doi errors out, just use preprint doi
            return None, None
