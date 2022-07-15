from backend.gene_info.config import GeneInfoConfig
import logging
import urllib.request
import json
import xml.etree.ElementTree as ET


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
        try:
            self.api_key = f"&api_key={gene_info_config.ncbi_api_key}"
        except RuntimeError:
            logging.error("Could not find NCBI API key")
            self.api_key = None

    def fetch_gene_info_tree(self, uid):
        """
        Given a gene UID from NCBI, returns an XML tree of gene information
        """

        if self.api_key is None:
            raise NCBIAPIException
        fetch_url = f"{self.base_ncbi_uri}efetch.fcgi?db=gene&id={uid}{self.api_key}&retmode=xml"
        try:
            return urllib.request.urlopen(fetch_url).read()
        except Exception:
            raise NCBIUnexpectedResultException

    def fetch_gene_uid(self, gene, geneID):
        """
        Given a gene ensembl ID, returns NCBI's corresponding gene UID
        """

        if self.api_key is None:
            raise NCBIAPIException
        search_url = f"{self.base_ncbi_uri}esearch.fcgi?db=gene&term={geneID}{self.api_key}&retmode=json"
        try:
            search_response = urllib.request.urlopen(search_url).read()
        except Exception:
            raise NCBIUnexpectedResultException
        search_result = self._load_search_result(search_response)
        if not search_result:
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
            logging.error(f"Unexpected NCBI search result, got {search_result}")
            return None
        return int(search_result["esearchresult"]["idlist"][0])

    def parse_gene_info_tree(self, tree_response):
        # parse NCBI XML response into relevant values to return by gene_info API
        result_tree = ET.ElementTree(ET.fromstring(tree_response))
        root = result_tree.getroot()
        synonyms = []
        summary = ""
        name = ""

        summary_tag = "Entrezgene_summary"
        gene_tag = "Entrezgene_gene"
        desc_tag = "Gene-ref_desc"
        syn_tag = "Gene-ref_syn"

        if len(root) > 0:
            for x in root[0]:
                if x.tag == summary_tag:
                    summary = x.text
                elif x.tag == gene_tag:
                    if len(x) > 0:
                        for y in x[0]:
                            if y.tag == desc_tag:
                                name = y.text
                            elif y.tag == syn_tag:
                                for syn in y:
                                    synonyms.append(syn.text)
        return dict(
            name=name,
            summary=summary,
            synonyms=synonyms,
        )
