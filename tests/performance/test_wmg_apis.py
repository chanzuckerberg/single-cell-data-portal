import json
import logging
import os
import timeit
import unittest

import requests

from tests.functional.backend.wmg.fixtures import (
    secondary_filter_common_case_request_data,
    secondary_filter_extreme_case_request_data,
)

# Note that these tests share fixtures and general test paths with the wmg api functional tests

logger = logging.getLogger(__name__)


@unittest.skipIf(os.getenv("DEPLOYMENT_STAGE") != "prod", "this test should only run in prod")
@unittest.skip("Skipping WMG V1 Performance Tests. WMG V1 API is deprecated. These tests will be ported to WMG V2")
class TestWmgApiPerformanceProd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.api = "https://api.cellxgene.cziscience.com"
        cls.api = f"{cls.api}/wmg/v1"

        # get snapshot id
        res = requests.get(f"{cls.api}/primary_filter_dimensions")
        cls.data = json.loads(res.content)

    def test_primary_filters(self):
        """
        Load primary filters in less than 1 second
        """
        MAX_RESPONSE_TIME_SECONDS = 1

        def make_request():
            return requests.get(f"{self.api}/primary_filter_dimensions")

        seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request, number=10)
        logger.info(f"primary filters seconds_for_10_runs: {seconds_for_10_runs}")
        self.assertGreater(MAX_RESPONSE_TIME_SECONDS * 10, seconds_for_10_runs)

    def test_secondary_filters_common_case(self):
        """
        1 tissue w/50 cell types, 20 genes, 3 secondary filters specified
        Returns in less than 10 seconds
        """
        MAX_RESPONSE_TIME_SECONDS = 10
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]

        def make_request():
            return requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)

        seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request, number=10)
        logger.info(f"secondary filters seconds_for_10_runs: {seconds_for_10_runs}")
        self.assertGreater(MAX_RESPONSE_TIME_SECONDS * 10, seconds_for_10_runs)

    def test_secondary_filters_extreme_case(self):
        """
        4 tissues w/largest cell type counts, 400 genes, no secondary filtering
        Returns in less than 15 seconds
        """
        MAX_RESPONSE_TIME_SECONDS = 15
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]

        def make_request():
            return requests.post(f"{self.api}/query", data=json.dumps(data), headers=headers)

        seconds_for_10_runs = timeit.timeit(setup="", stmt=make_request, number=10)
        logger.info(f"secondary filters extreme case: seconds_for_10_runs: {seconds_for_10_runs}")
        self.assertGreater(MAX_RESPONSE_TIME_SECONDS * 10, seconds_for_10_runs)
