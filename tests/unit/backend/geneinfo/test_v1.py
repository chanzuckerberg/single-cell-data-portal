import unittest
import requests
import json

from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from backend.corpora.api_server.app import app
from tests.unit.backend.geneinfo.fixtures import (
    correct_ensembl_to_gene_mappings
)


class GeneInfoAPIv1Tests(unittest.TestCase):
    ''' Tests Gene Info API endpoint '''

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    def test_correct_geneids(self):
        """
        Successfully returns correct response for ensembl IDs
        Includes responses that are missing information (returns empty strings)
        """
        for id in correct_ensembl_to_gene_mappings.keys():
            res = self.app.get(f"/geneinfo/v1/geneinfo?geneID={id}")
            data = json.loads(res.data)
            self.assertEqual(res.status_code, requests.codes.ok)
            self.assertEqual(data, correct_ensembl_to_gene_mappings[id])

    def test_incorrect_geneids(self):
        """
        Successfully returns 404 not found for ensembl IDs that do not exist
        """
        res1 = self.app.get("/geneinfo/v1/geneinfo?geneID=")
        data1 = json.loads(res1.data)
        self.assertEqual(res1.status_code, requests.codes.not_found)
        self.assertEqual(data1, "Unexpected NCBI search result")

        res2 = self.app.get("/geneinfo/v1/geneinfo?geneID=abc")
        data2 = json.loads(res2.data)
        self.assertEqual(res2.status_code, requests.codes.not_found)
        self.assertEqual(data2, "Unexpected NCBI search result")
