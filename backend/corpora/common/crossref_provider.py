import requests
from .corpora_config import CorporaConfig


class CrossrefProvider(object):
    """
    Provider class used to call Crossref and retrieve publisher metadata
    """

    def __init__(self, crossref_api_uri = None) -> None:
        # crossref_api_uri can be passed to the constructor and overrides the configuration for testing
        try:
            self.base_crossref_uri = crossref_api_uri or CorporaConfig().crossref_api_uri
        except RuntimeError:
            self.base_crossref_uri = None
        super().__init__()

    @staticmethod
    def parse_date_parts(obj):
        date_parts = obj["date-parts"][0]
        return (date_parts[0], date_parts[1])

    def fetch_metadata(self, doi):
        """
        Fetches and extracts publisher metadata from Crossref for a specified DOI.
        If the Crossref API URI isn't in the configuration, we will just return an empty object.
        This is to avoid calling Crossref in non-production environments.
        """
        if self.base_crossref_uri is not None:
            # TODO: if we're using the commercial API, the token should also be parametrized
            res = requests.get(f"{self.base_crossref_uri}/{doi}")
            message = res["message"]

            # Date
            published_date = (
                message.get("published-print") or message.get("published") or message.get("published-online")
            )
            published_year, published_month = self.parse_date_parts(published_date)

            # Journal
            if "short-container-title" in message and message["short-container-title"]:
                journal = message["short-container-title"][0]
            elif "container-title" in message and message["container-title"]:
                journal = message["container-title"][0]
            elif "institution" in message:
                journal = message["institution"][0]["name"]
            else:
                journal = None

            # Authors
            authors = message["author"]
            parsed_authors = [{"given": a["given"], "family": a["family"]} for a in authors]

            # Preprint
            is_preprint = message.get("subtype") == "preprint"

            return {
                "authors": parsed_authors,
                "published_year": published_year,
                "published_month": published_month,
                "journal": journal,
                "is_preprint": is_preprint,
            }

        else:
            return None
