import json
import logging
import urllib.request
import xml.etree.ElementTree as ET

from backend.gene_info.config import GeneInfoConfig


class NCBIException(Exception):
    pass


class NCBIAPIException(NCBIException):
    pass


class NCBIUnexpectedResultException(NCBIException):
    pass


class NCBIProvider:
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
            self.api_key = None  # type: ignore

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
            raise NCBIUnexpectedResultException from None

    def fetch_gene_uid(self, geneID, gene):
        """
        Given a gene ensembl ID and gene name, returns a tuple with NCBI's corresponding gene UID and a
        boolean noting if the result is from searching by the gene name, instead of the gene ENSEMBL id.
        Initially, uses ensembl ID to find gene UID, but in the event that this returns an
        unexpected result, call the NCBI API again with gene name. This is successful if
        the response returns only 1 result for UID.
        """

        if self.api_key is None:
            raise NCBIAPIException
        # first search with gene ENSEMBL id
        try:
            search_response = self._search_gene_uid(geneID)
        except NCBIUnexpectedResultException:
            raise NCBIUnexpectedResultException from None

        # search with gene name if needed
        if not self._is_valid_search_result(search_response) and gene and gene != "":
            try:
                # to refine the search parameters, gene is searched with the label "Gene Name" and human "Organism"
                gene_search_term = "(" + str(gene) + "%5BGene%20Name%5D)%20AND%20human%5BOrganism%5D"
                search_response = self._search_gene_uid(gene_search_term)
                if self._is_valid_search_result(search_response):
                    show_warning_banner = False
                    return (int(search_response["esearchresult"]["idlist"][0]), show_warning_banner)
                else:
                    logging.error(f"Unexpected NCBI search result, got {search_response}")
                    raise NCBIUnexpectedResultException
            except NCBIUnexpectedResultException:
                raise NCBIUnexpectedResultException from None
        elif not self._is_valid_search_result(search_response):
            logging.error(f"Unexpected NCBI search result, got {search_response}")
            raise NCBIUnexpectedResultException
        else:
            show_warning_banner = False
            return (int(search_response["esearchresult"]["idlist"][0]), show_warning_banner)

    def _search_gene_uid(self, term):
        """
        Conducts an Esearch using NCBI's E-Utilities API with provided term
        """
        search_url = f"{self.base_ncbi_uri}esearch.fcgi?db=gene&term={term}{self.api_key}&retmode=json"
        try:
            search_response = urllib.request.urlopen(search_url).read()
        except Exception:
            raise NCBIUnexpectedResultException from None
        return json.loads(search_response)

    def _is_valid_search_result(self, search_result):
        """
        Checks that a search result contains only one UID as a result
        """
        try:
            int(search_result["esearchresult"]["idlist"][0])
            if len(search_result["esearchresult"]["idlist"]) != 1:
                return False
        except (ValueError, KeyError, IndexError):
            return False
        return True

    def parse_gene_info_tree(self, tree_response):
        """
        Parse NCBI XML response into relevant values to return by gene_info API
        """
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
                elif x.tag == gene_tag and len(x) > 0:
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
