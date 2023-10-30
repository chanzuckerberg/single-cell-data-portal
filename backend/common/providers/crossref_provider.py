import logging
from datetime import datetime
from urllib.parse import urlparse

import requests

from backend.common.corpora_config import CorporaConfig


class CrossrefException(Exception):
    pass


class CrossrefFetchException(CrossrefException):
    pass


class CrossrefDOINotFoundException(CrossrefException):
    pass


class CrossrefParseException(CrossrefException):
    pass


class CrossrefProvider:
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
        # Remove the https://doi.org part
        parsed = urlparse(doi)
        if parsed.scheme and parsed.netloc:
            doi = parsed.path

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

        return self.fetch_preprint_published_doi(res)

    def fetch_metadata(self, doi: str) -> dict:
        """
        Fetches and extracts publisher metadata from Crossref for a specified DOI.
        If the Crossref API URI isn't in the configuration, we will just return an empty object.
        This is to avoid calling Crossref in non-production environments.
        """

        res = self._fetch_crossref_payload(doi)
        if not res:
            return

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
                    journal = message["short-container-title"][0]
                elif "container-title" in message and message["container-title"]:
                    journal = message["container-title"][0]
                elif "institution" in message:
                    journal = message["institution"][0]["name"]
            except Exception:
                raise CrossrefParseException("Journal node missing") from None

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

            return {
                "authors": parsed_authors,
                "published_year": published_year,
                "published_month": published_month,
                "published_day": published_day,
                "published_at": datetime.timestamp(datetime(published_year, published_month, published_day)),
                "journal": journal,
                "is_preprint": is_preprint,
            }
        except Exception as e:
            raise CrossrefParseException("Cannot parse metadata from Crossref") from e

    def fetch_preprint_published_doi(self, res):
        """
        If CrossRef API Response is of a preprint with a known published DOI, fetch the metadata for the published DOI
        instead. Otherwise, return the input Response.
        """
        message = res.json()["message"]
        is_preprint = message.get("subtype") == "preprint"
        has_published_doi = ("relation" in message) and ("is-preprint-of" in message["relation"])

        if is_preprint and has_published_doi:
            published_ids = message["relation"]["is-preprint-of"]
            for pub_id in published_ids:
                if pub_id["id-type"] == "doi":
                    published_doi = published_ids[pub_id]["id"].replace("\\", "")
                    try:
                        published_doi_res = requests.get(
                            f"{self.base_crossref_uri}/{published_doi}",
                            headers={"Crossref-Plus-API-Token": f"Bearer {self.crossref_api_key}"},
                        )
                        published_doi_res.raise_for_status()
                        return published_doi_res
                    except requests.RequestException:
                        pass

        # return pre-print if fetching original doi fails
        return res
