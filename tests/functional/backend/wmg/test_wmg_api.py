import json
import unittest

import requests  # type: ignore

from tests.functional.backend.common import BaseFunctionalTestCase
from tests.functional.backend.wmg.fixtures import (
    markers_happy_path,
    markers_missing_tissue,
    secondary_filter_common_case_request_data,
    secondary_filter_data_with_ontology_term_ids,
    secondary_filter_extreme_case_request_data,
)

# Note that these tests share fixtures and general test paths with the wmg api performance tests


@unittest.skip("Skipping WMG V1 Functional Tests. WMG V1 API is deprecated. These tests will be ported to WMG V2")
class TestWmgApi(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.api = f"{cls.api}/wmg/v1"

        # get snapshot id
        res = requests.get(f"{cls.api}/primary_filter_dimensions")
        cls.data = json.loads(res.content)

    def test_primary_filters(self):
        """
        Load primary filters in less than 1.5 seconds
        """
        res = self.session.get(f"{self.api}/primary_filter_dimensions")
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_query_endpoint_common_case(self):
        """
        1 tissue w/50 cell types, 20 genes, 3 secondary filters specified
        Returns in less than 10 seconds
        """
        headers = {"Content-Type": "application/json"}

        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/query", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_query_endpoint_extreme_case(self):
        """
        4 tissues w/largest cell type counts, 400 genes, no secondary filtering
        Returns in less than 15 seconds
        """
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/query", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_filter_endpoint_common_case(self):
        """
        /v1/filters should support the common case /v1/queries supports
        """
        headers = {"Content-Type": "application/json"}

        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/filters", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_filter_endpoint_extreme_case(self):
        """
        /v1/filters should support the extreme case /v1/queries supports
        """
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/filters", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_filter_endpoint_supports_ontology_term_ids(self):
        """
        /v1/filters differs from /v1/query in that it supports the cell_type_ontology_term_ids filter
        Ensure that hitting this endpoint with cell_type_ontology_term_ids is a valid request
        """
        headers = {"Content-Type": "application/json"}

        data = secondary_filter_data_with_ontology_term_ids.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/filters", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_markers_happy_path(self):
        headers = {"Content-Type": "application/json"}

        data = markers_happy_path.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/markers", data=json.dumps(data), headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertGreater(len(res.content), 10)

    def test_markers_missing_tissue(self):
        """
        requests missing required params should fail
        """
        headers = {"Content-Type": "application/json"}

        data = markers_missing_tissue.copy()
        data["snapshot_id"] = self.data["snapshot_id"]
        res = self.session.post(f"{self.api}/markers", data=json.dumps(data), headers=headers)
        self.assertStatusCode(400, res)
        self.assertGreater(len(res.content), 10)
