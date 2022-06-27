from ..config import GeneInfoConfig
import logging
import urllib.request
import json


class NCBIException(Exception):
    pass


class NCBIAPIException(NCBIException):
    pass


class NCBIUnexpectedResultException(NCBIException):
    pass


class NCBIProvider(object):
    """
    Provider class used to generate NCBI URL
    """

    def __init__(self) -> None:
        self.base_ncbi_uri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        gene_info_config = GeneInfoConfig()
        gene_info_config.load()
        try:
            self.ncbi_api_key = f"&api_key={gene_info_config.ncbi_api_key}"
        except RuntimeError:
            logging.error("Could not find NCBI API key")
            self.ncbi_api_key = None
        super().__init__()

    def _get_api_key(self):
        if self.ncbi_api_key is None:
            return None
        else:
            return f"&api_key={self.ncbi_api_key}"

    def fetch_gene_info_tree(self, uid):
        """
        Given a gene UID from NCBI, returns an XML tree of gene information
        """

        if self._get_api_key() is None:
            raise NCBIAPIException
        fetch_url = f"{self.base_ncbi_uri}efetch.fcgi?db=gene&id={uid}{self._get_api_key()}&retmode=xml"
        try:
            return urllib.request.urlopen(fetch_url).read()
        except Exception:
            raise NCBIUnexpectedResultException

    def fetch_gene_uid(self, geneID):
        """
        Given a gene ensembl ID, returns NCBI's corresponding gene UID
        """

        if self._get_api_key() is None:
            raise NCBIAPIException
        search_url = f"{self.base_ncbi_uri}esearch.fcgi?db=gene&term={geneID}{self._get_api_key()}&retmode=json"
        try:
            search_response = urllib.request.urlopen(search_url).read()
        except Exception:
            raise NCBIUnexpectedResultException
        search_result = self._load_search_result(search_response)
        if not search_result:
            logging.error(f"Unexpected NCBI search result, got {search_result}")
            raise NCBIUnexpectedResultException
        else:
            return search_result

    def _load_search_result(self, response):
        search_result = json.loads(response)
        if (
            "esearchresult" not in search_result
            or "idlist" not in search_result["esearchresult"]
            or len(search_result["esearchresult"]["idlist"]) < 1
            or int(search_result["esearchresult"]["count"]) != 1
        ):
            return None
        return int(search_result["esearchresult"]["idlist"][0])
