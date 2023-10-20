import json
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import call, patch

import requests

from backend.api_server.app import app
from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.gene_info.api import ensembl_ids, ncbi_provider
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup


class GeneInfoAPIv1Tests(unittest.TestCase):
    """Tests Gene Info API endpoint"""

    def setUp(self):
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)
        self.final_gene_info_result = {
            "name": "",
            "summary": "",
            "ncbi_url": "https://www.ncbi.nlm.nih.gov/gene/1",
            "synonyms": [],
            "show_warning_banner": False,
        }

        self.mock_config = patch("backend.gene_info.config.GeneInfoConfig.__getattr__", return_value="mock_key")
        self.mock_config.start()

    @classmethod
    def setUpClass(cls) -> None:
        cls.maxDiff = None

    def tearDown(self):
        self.mock_config.stop()

    @patch("backend.gene_info.api.ncbi_provider.urllib.request.urlopen")
    @patch("backend.gene_info.api.v1.NCBIProvider._search_gene_uid")
    def test_api_fetch_calls(self, mock_search_gene_uid, mock_get):
        """
        Mocks API key for NCBI requests, checks call counts and the correct external API calls.
        Will break if the external API call is down!
        """
        mock_get.read = None
        mock_search_gene_uid.return_value = {"esearchresult": {"idlist": [348]}}
        test_fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=348&retmode=xml"
        test_provider = ncbi_provider.NCBIProvider()
        test_provider.api_key = ""
        test_provider.fetch_gene_uid("ENSG00000130203", "APOE")
        test_provider.fetch_gene_info_tree(348)
        mock_get.assert_has_calls([call(test_fetch_url), call().read()])
        self.assertEqual(mock_get.call_count, 1)
        self.assertEqual(mock_search_gene_uid.call_count, 1)

    @patch("backend.gene_info.api.ncbi_provider.urllib.request.urlopen")
    def test_api_search_calls(self, mock_get):
        """
        In combination with test_api_fetch_calls, checks external API call in search requests (_search_gene_uid)
        """

        def mock_json_loads(value):
            return value

        with patch("backend.gene_info.api.ncbi_provider.json.loads", mock_json_loads):
            test_provider = ncbi_provider.NCBIProvider()
            test_provider.api_key = ""
            test_search_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                "db=gene&term=ENSG00000130203&retmode=json"
            )
            test_provider._search_gene_uid("ENSG00000130203")
            mock_get.assert_has_calls([call(test_search_url), call().read()])
            self.assertEqual(mock_get.call_count, 1)

    @patch("backend.gene_info.api.v1.NCBIProvider._search_gene_uid")
    @patch("backend.gene_info.api.v1.NCBIProvider._is_valid_search_result")
    def test_api_double_calls(self, mock_is_valid_search_result, mock_search_gene_uid):
        """Checks that, if a gene fails on the first ID search, ncbi_provider will search again (with gene name)"""
        mock_search_gene_uid.return_value = "failed"
        mock_is_valid_search_result.return_value = False
        self.app.get("/gene_info/v1/gene_info?geneID=ensembl&gene=name")
        self.assertEqual(mock_search_gene_uid.call_count, 2)
        mock_search_gene_uid.assert_has_calls(
            [call("ensembl"), call("(name%5BGene%20Name%5D)%20AND%20human%5BOrganism%5D")]
        )
        self.assertEqual(mock_is_valid_search_result.call_count, 2)
        mock_is_valid_search_result.assert_has_calls([call("failed"), call("failed")])
        self.assertRaises(NotFoundHTTPException)

    def test_unexpected_ncbi_responses(self):
        """Correctly identifies a valid search result from NCBI"""
        test_provider = ncbi_provider.NCBIProvider()
        self.assertFalse(test_provider._is_valid_search_result({}))
        self.assertFalse(test_provider._is_valid_search_result({"esearchresult": {}}))
        self.assertFalse(test_provider._is_valid_search_result({"esearchresult": {"idlist": []}}))
        self.assertFalse(test_provider._is_valid_search_result({"esearchresult": {"idlist": ["hello"]}}))
        self.assertFalse(test_provider._is_valid_search_result({"esearchresult": {"idlist": ["1", "2"]}}))
        self.assertTrue(test_provider._is_valid_search_result({"esearchresult": {"idlist": ["1"]}}))

    @patch("backend.gene_info.api.v1.NCBIProvider.fetch_gene_uid")
    @patch("backend.gene_info.api.v1.NCBIProvider.fetch_gene_info_tree")
    @patch("backend.gene_info.api.v1.NCBIProvider.parse_gene_info_tree")
    @patch("backend.gene_info.api.v1.GeneChecker.get_id")
    def test_gene_info_with_id(
        self, mock_get_id, mock_parse_gene_info_tree, mock_fetch_gene_info_tree, mock_fetch_gene_uid
    ):
        """
        Successfully calls NCBIProvider fetch and search functions with a parameter of only gene ID
        """
        mock_fetch_gene_uid.return_value = (1, False)
        mock_fetch_gene_info_tree.return_value = None
        mock_parse_gene_info_tree.return_value = self.final_gene_info_result
        mock_get_id.return_value = "ensembl1"

        res = self.app.get("/gene_info/v1/gene_info?geneID=ensembl1")
        self.assertTrue(mock_fetch_gene_uid.called)
        self.assertTrue(mock_fetch_gene_info_tree.called)
        self.assertTrue(mock_parse_gene_info_tree.called)
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertEqual(json.loads(res.data), self.final_gene_info_result)

    @patch("backend.gene_info.api.v1.NCBIProvider.fetch_gene_uid")
    @patch("backend.gene_info.api.v1.NCBIProvider.fetch_gene_info_tree")
    @patch("backend.gene_info.api.v1.NCBIProvider.parse_gene_info_tree")
    @patch("backend.gene_info.api.v1.GeneChecker.get_id")
    def test_gene_info_with_id_and_name(
        self, mock_get_id, mock_parse_gene_info_tree, mock_fetch_gene_info_tree, mock_fetch_gene_uid
    ):
        """
        Successfully calls NCBIProvider fetch and search functions with both ID and name params
        """
        mock_fetch_gene_uid.return_value = (1, False)
        mock_fetch_gene_info_tree.return_value = None
        mock_parse_gene_info_tree.return_value = self.final_gene_info_result
        mock_get_id.return_value = "ensembl1"

        res = self.app.get("/gene_info/v1/gene_info?geneID=ensembl1&gene=name")
        self.assertTrue(mock_fetch_gene_uid.called)
        self.assertTrue(mock_fetch_gene_info_tree.called)
        self.assertTrue(mock_parse_gene_info_tree.called)
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertEqual(json.loads(res.data), self.final_gene_info_result)

    @patch("backend.gene_info.api.v1.GeneChecker.get_id")
    def test_gene_checker_call(self, mock_get_id):
        """
        Successfully calls GeneChecker if called with only a name parameter
        """
        self.app.get("/gene_info/v1/gene_info?gene=APOE")
        self.assertTrue(mock_get_id.called)

    @patch("backend.gene_info.api.v1.gene_info")
    def test_incorrect_gene_ids(self, test_provider):
        """
        Successfully raises exception for ensembl IDs that do not exist
        """
        test_provider.provider.api_key = ""
        res1 = self.app.get("/gene_info/v1/gene_info?geneID=")
        self.assertEqual(res1.status_code, 404)
        self.assertEqual(json.loads(res1.data)["detail"], "No gene label provided.")
        res2 = self.app.get("/gene_info/v1/gene_info?geneID=abc")
        self.assertEqual(res2.status_code, 404)
        self.assertEqual(json.loads(res2.data)["detail"], "Unexpected NCBI search result")

    def test_correct_parse_xml_tree(self):
        """
        Successfully parses a properly formatted XML tree for gene information
        """
        provider = ncbi_provider.NCBIProvider()
        root = ET.Element("root")
        entrez = ET.SubElement(root, "Entrez")
        summary = ET.SubElement(entrez, "Entrezgene_summary")
        summary.text = "gene summary"
        gene_root = ET.SubElement(entrez, "Entrezgene_gene")
        gene = ET.SubElement(gene_root, "gene_root")
        desc = ET.SubElement(gene, "Gene-ref_desc")
        desc.text = "gene description"
        syn = ET.SubElement(gene, "Gene-ref_syn")
        syn1 = ET.Element("syn1")
        syn1.text = "syn1"
        syn.append(syn1)
        syn2 = ET.Element("syn2")
        syn2.text = "syn2"
        syn.append(syn2)
        self.assertEqual(
            provider.parse_gene_info_tree(ET.tostring(root)),
            dict(name="gene description", summary="gene summary", synonyms=["syn1", "syn2"]),
        )

    def test_parse_xml_tree_missing_values(self):
        """
        Returns correct dictionary with missing values if given XML tree with missing values
        """
        provider = ncbi_provider.NCBIProvider()

        root1 = ET.Element("root")
        self.assertEqual(provider.parse_gene_info_tree(ET.tostring(root1)), dict(name="", summary="", synonyms=[]))

        root2 = ET.Element("root")
        ET.SubElement(root2, "Entrez")
        self.assertEqual(provider.parse_gene_info_tree(ET.tostring(root2)), dict(name="", summary="", synonyms=[]))

        root3 = ET.Element("root")
        entrez3 = ET.SubElement(root3, "Entrez")
        summary3 = ET.SubElement(entrez3, "Entrezgene_summary")
        summary3.text = "gene summary"
        self.assertEqual(
            provider.parse_gene_info_tree(ET.tostring(root3)), dict(name="", summary="gene summary", synonyms=[])
        )

    def test_gene_checker(self):
        """
        GeneChecker successfully creates a dictionary of gene names to gene IDs,
        checks validity of gene names, and returns ensembl ID for a gene name.
        """
        gene_checker = ensembl_ids.GeneChecker()
        self.assertTrue(gene_checker.gene_dict)
        self.assertTrue("APOE" in gene_checker.gene_dict)
        self.assertFalse("not a gene" in gene_checker.gene_dict)
        self.assertEqual(gene_checker.get_id("APOE"), "ENSG00000130203")
