import html
import logging
from datetime import datetime
from typing import Optional, Tuple

import requests

from backend.common.citation import format_citation_crossref
from backend.common.corpora_config import CorporaConfig
from backend.common.doi import doi_curie_from_link


class CrossrefProviderInterface:
    def fetch_metadata(self, doi: str) -> Tuple[Optional[dict], Optional[str], Optional[float]]:
        return None, None, None

    def fetch_preprint_published_doi(self, doi):
        pass

    def _fetch_crossref_payload(self, doi):
        pass

    def get_title_and_citation_from_doi(self, doi: str) -> str:
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

    def get_title_and_citation_from_doi(self, doi: str) -> str:
        """
        Retrieves the title and citation from a DOI.

        Parameters
        ----------
        doi : str
            The DOI string.

        Returns
        -------
        str
            The title and citation associated with the DOI.
        """

        response = self._fetch_crossref_payload(doi)
        data = response.json()

        # Get the title and citation count from the data
        try:
            title = data["message"]["title"][0]
            citation = format_citation_crossref(data["message"])
        except Exception:
            try:
                title = data["message"]["items"][0]["title"][0]
                citation = format_citation_crossref(data["message"]["items"][0])
            except Exception:
                return doi
        return f"{title}\n\n - {citation}"

    def fetch_metadata(self, doi: str) -> Tuple[Optional[dict], Optional[str], Optional[datetime]]:
        """
        Fetches and extracts publisher metadata from Crossref for a specified DOI.
        If the Crossref API URI isn't in the configuration, we will just return an empty object.
        This is to avoid calling Crossref in non-production environments.
        :param doi: str - DOI uri link or curie identifier
        return: tuple - publisher metadata dict and DOI curie identifier
        """

        doi_curie = doi_curie_from_link(doi)

        res = self._fetch_crossref_payload(doi_curie)
        if not res:
            return None, None, None

        try:
            message = res.json()["message"]

            # Date
            published_date = (
                message.get("published-print") or message.get("published") or message.get("published-online")
            )

            if published_date is None:
                raise CrossrefParseException("Date node missing")

            published_year, published_month, published_day = self.parse_date_parts(published_date)

            # Calculate the deposited date; used when checking for updates.
            deposited_at = None
            if "deposited" in message and (deposited_timestamp := message["deposited"].get("timestamp")) is not None:
                deposited_at = deposited_timestamp / 1000

            # Journal
            try:
                raw_journal = None
                if "short-container-title" in message and message["short-container-title"]:
                    raw_journal = message["short-container-title"][0]
                elif "container-title" in message and message["container-title"]:
                    raw_journal = message["container-title"][0]
                elif "institution" in message:
                    raw_journal = message["institution"][0]["name"]

                if raw_journal is None:
                    raise CrossrefParseException("Journal node missing")
            except Exception:
                raise CrossrefParseException("Journal node missing") from None

            journal = html.unescape(raw_journal)

            # Authors
            # Note: make sure that the order is preserved, as it is a relevant information
            authors = message["author"]
            parsed_authors = []
            for author in authors:
                if "given" in author and "family" in author:
                    parsed_author = {"given": author["given"], "family": author["family"]}
                    if parsed_author not in parsed_authors:
                        parsed_authors.append(parsed_author)
                elif "family" in author:
                    # Assume family is consortium
                    parsed_authors.append({"name": author["family"]})
                elif "name" in author:
                    parsed_authors.append({"name": author["name"]})

            # Preprint
            is_preprint = message.get("subtype") == "preprint"
            if is_preprint:
                published_metadata, published_doi_curie, published_deposited_at = self.fetch_published_metadata(message)
                if published_metadata and published_doi_curie:  # if not, use preprint doi curie
                    return published_metadata, published_doi_curie, published_deposited_at

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
                deposited_at,
            )
        except Exception as e:
            raise CrossrefParseException("Cannot parse metadata from Crossref") from e

    def fetch_published_metadata(
        self, doi_response_message: dict
    ) -> Tuple[Optional[dict], Optional[str], Optional[float]]:
        try:
            published_doi = doi_response_message["relation"]["is-preprint-of"]
            # the new DOI to query for ...
            for entity in published_doi:
                if entity["id-type"] == "doi":
                    return self.fetch_metadata(entity["id"])
        except Exception:  # if fetch of published doi errors out, just use preprint doi
            return None, None, None
