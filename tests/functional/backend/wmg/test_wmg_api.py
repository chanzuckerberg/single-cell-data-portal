import json

import pytest
import requests

from tests.functional.backend.wmg.fixtures import (
    markers_happy_path,
    markers_missing_tissue,
    secondary_filter_common_case_request_data,
    secondary_filter_data_with_ontology_term_ids,
    secondary_filter_extreme_case_request_data,
)


@pytest.fixture(scope="module")
def setup(api_url):
    api = f"{api_url}/wmg/v1"
    res = requests.get(f"{api}/primary_filter_dimensions")
    primary_filter_dimensions = json.loads(res.content)
    return api, primary_filter_dimensions


# Note that these tests share fixtures and general test paths with the wmg api performance tests


@pytest.mark.skip("Skipping WMG V1 Functional Tests. WMG V1 API is deprecated. These tests will be ported to WMG V2")
class TestWmgApi:
    def test_primary_filters(self, setup, session):
        """
        Load primary filters in less than 1.5 seconds
        """
        api, _ = setup
        res = session.get(f"{api}/primary_filter_dimensions")
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_query_endpoint_common_case(self, setup, session):
        """
        1 tissue w/50 cell types, 20 genes, 3 secondary filters specified
        Returns in less than 10 seconds
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/query", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_query_endpoint_extreme_case(self, setup, session):
        """
        4 tissues w/largest cell type counts, 400 genes, no secondary filtering
        Returns in less than 15 seconds
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/query", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_filter_endpoint_common_case(self, setup, session):
        """
        /v1/filters should support the common case /v1/queries supports
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_common_case_request_data.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/filters", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_filter_endpoint_extreme_case(self, setup, session):
        """
        /v1/filters should support the extreme case /v1/queries supports
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}
        data = secondary_filter_extreme_case_request_data.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{self.api}/filters", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_filter_endpoint_supports_ontology_term_ids(self, setup, session):
        """
        /v1/filters differs from /v1/query in that it supports the cell_type_ontology_term_ids filter
        Ensure that hitting this endpoint with cell_type_ontology_term_ids is a valid request
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}

        data = secondary_filter_data_with_ontology_term_ids.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/filters", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_markers_happy_path(self, setup, session):
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}

        data = markers_happy_path.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/markers", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.ok
        assert len(res.content) > 10

    def test_markers_missing_tissue(self, setup, session):
        """
        requests missing required params should fail
        """
        api, primary_filter_dimensions = setup
        headers = {"Content-Type": "application/json"}

        data = markers_missing_tissue.copy()
        data["snapshot_id"] = primary_filter_dimensions["snapshot_id"]
        res = session.post(f"{api}/markers", data=json.dumps(data), headers=headers)
        assert res.status_code == requests.codes.bad_request
        assert len(res.content) > 10
