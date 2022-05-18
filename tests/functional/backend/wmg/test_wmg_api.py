import json
import os
import timeit
import unittest

import requests

from tests.functional.backend.common import API_URL
from tests.functional.backend.wmg.fixtures import genes_20_count, genes_400_count


class TestWmgApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]
        cls.api = API_URL.get(cls.deployment_stage)
        cls.api = f"{cls.api}/wmg/v1"

        # get snapshot id
        res = requests.get(f"{cls.api}/primary_filter_dimensions")
        cls.data = json.loads(res.content)

    def test_primary_filters(self):
        """
        Load primary filters in less than 1.5 seconds
        """
        MAX_RESPONSE_TIME_SECONDS = 1.5
        with self.subTest("Run once"):
            res = requests.get(f"{self.api}/primary_filter_dimensions")
            self.assertGreater(MAX_RESPONSE_TIME_SECONDS, res.elapsed.seconds)
            self.assertEqual(res.status_code, requests.codes.ok)
            self.assertGreater(len(res.content), 10)
        with self.subTest("Run 10 times"):
            def make_request():
                return requests.get(f"{self.api}/primary_filter_dimensions")

            seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request, number=10)
            self.assertGreater(MAX_RESPONSE_TIME_SECONDS*10, seconds_for_10_runs)

    def test_secondary_filters_common_case(self):
        """
        1 tissue w/50 cell types, 20 genes, 3 secondary filters specified
        Returns in less than 10 seconds
        """
        MAX_RESPONSE_TIME_SECONDS = 10
        headers = {"Content-Type": "application/json"}
        data = {
            "filter": {
                "dataset_ids": [],
                "disease_ontology_term_ids": ["MONDO:0005812", "MONDO:0100096", "PATO:0000461"],
                "ethnicity_ontology_term_ids": ["HANCESTRO:0010", "HANCESTRO:0014"],
                "gene_ontology_term_ids": list(genes_20_count.keys()),
                "organism_ontology_term_id": "NCBITaxon:9606",
                "sex_ontology_term_ids": ["PATO:0000383"],
                "tissue_ontology_term_ids": ["UBERON:0000178"],
            },  # blood (more than 50 cell types)
            "include_filter_dims": True,
            "snapshot_id": self.data["snapshot_id"],
        }
        with self.subTest("Run once"):
            res = requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)
            self.assertGreaterEqual(MAX_RESPONSE_TIME_SECONDS, res.elapsed.seconds)
            self.assertEqual(res.status_code, requests.codes.ok)
            self.assertGreater(len(res.content), 10)
        with self.subTest("Run 10 times"):
            with self.subTest("Test 10 times"):
                def make_request():
                    return requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)

                seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request, number=10)
                self.assertGreater(MAX_RESPONSE_TIME_SECONDS * 10, seconds_for_10_runs)

    def test_secondary_filters_extreme_case(self):
        """
        4 tissues w/largest cell type counts, 400 genes, no secondary filtering
        Returns in less than 15 seconds
        """
        MAX_RESPONSE_TIME_SECONDS = 15
        headers = {"Content-Type": "application/json"}
        data = {
            "filter": {
                "dataset_ids": [],
                "disease_ontology_term_ids": [],
                "ethnicity_ontology_term_ids": [],
                "gene_ontology_term_ids": list(genes_400_count.keys()),
                "organism_ontology_term_id": "NCBITaxon:9606",
                "sex_ontology_term_ids": [],
                "tissue_ontology_term_ids": [
                    "UBERON:0000178",
                    "UBERON:0002048",
                    "UBERON:0000029",
                    "UBERON:0000970",
                    "UBERON:0000362",
                ],
            },
            "include_filter_dims": True,
            "snapshot_id": self.data["snapshot_id"],
        }
        with self.subTest("Test once"):
            res = requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)
            self.assertGreater(MAX_RESPONSE_TIME_SECONDS, res.elapsed.seconds)
            self.assertEqual(res.status_code, requests.codes.ok)
            self.assertGreater(len(res.content), 10)

        with self.subTest("Test 10 times"):
            def make_request():
                return requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)

            seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request,
                                                number=10)
            self.assertGreater(MAX_RESPONSE_TIME_SECONDS*10, seconds_for_10_runs)
