import json
import os
import unittest

import requests

from tests.functional.backend.common import API_URL
from tests.functional.backend.wmg.fixtures import (
    secondary_filter_common_case_request_data,
    secondary_filter_extreme_case_request_data,
)


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
        res = requests.get(f"{self.api}/primary_filter_dimensions")
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertGreater(len(res.content), 10)

    def test_secondary_filters_common_case(self):
        """
        1 tissue w/50 cell types, 20 genes, 3 secondary filters specified
        Returns in less than 10 seconds
        """
        headers = {"Content-Type": "application/json"}

        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = requests.post(
            f"{self.api}/query", data=json.dumps(data), headers=headers
        )
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertGreater(len(res.content), 10)

    def test_secondary_filters_extreme_case(self):
        """
        4 tissues w/largest cell type counts, 400 genes, no secondary filtering
        Returns in less than 15 seconds
        """
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = requests.post(
            f"{self.api}/query", data=json.dumps(data), headers=headers
        )
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertGreater(len(res.content), 10)
