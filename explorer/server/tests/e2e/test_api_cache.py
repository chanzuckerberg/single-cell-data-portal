import json
import os
import unittest
from urllib.parse import quote

import requests


@unittest.skipIf(os.getenv("CXG_API_BASE") is None, "CXG_API_BASE not defined")
class TestAPICache(unittest.TestCase):
    """These tests must run against a deployed API backed by AWS cloudfront."""

    def setUp(self):
        self.api_url_base = os.environ["CXG_API_BASE"].rstrip("/")
        response = requests.get("/".join([self.api_url_base, "dp/v1/collections"]))
        response.raise_for_status()
        collection_id = response.json()["collections"][0]["id"]
        response = requests.get("/".join([self.api_url_base, f"dp/v1/collections/{collection_id}"]))
        response.raise_for_status()
        self.dataset_id = response.json()["datasets"][0]["id"]
        self.explorer_api_url_base = "/".join([self.api_url_base, "cellxgene", "e", f"{self.dataset_id}.cxg"])
        self.cloudfront_miss = "Miss from cloudfront"
        self.cloudfront_hit = "Hit from cloudfront"

    @staticmethod
    def generate_error_msg(response, message):
        return json.dumps(
            {"msg": message, "x-cache": response.headers["x-cache"], "x-amz-cf-id": response.headers["x-amz-cf-id"]}
        )

    def test_uncached_endpoints(self):
        """Verify the cache headers are properly set for uncached endpoints"""
        error_message = "Miss Expected"
        endpoints = ["api/v0.3/dataset-metadata", "api/v0.3/s3_uri"]

        def verify_cache_control(_response):
            _response.raise_for_status()
            cache_control = _response.headers["cache-control"]
            self.assertIn("no-store", cache_control)
            self.assertIn("max-age=0", cache_control)
            self.assertIn("public", cache_control)

        for endpoint in endpoints:
            with self.subTest(endpoint):
                # Test
                url = "/".join([self.explorer_api_url_base, endpoint])
                response = requests.head(url)

                # Verify
                verify_cache_control(response)
                self.assertEqual(
                    response.headers["x-cache"],
                    self.cloudfront_miss,
                    msg=self.generate_error_msg(response, error_message),
                )

                # Call the request twice to make sure the cache wasn't cold.
                response = requests.head(url)
                verify_cache_control(response)
                self.assertEqual(
                    response.headers["x-cache"],
                    self.cloudfront_miss,
                    msg=self.generate_error_msg(response, error_message),
                )

    def test_cached_endpoints(self):
        """Verify cloudfront caching headers are properly set for cached endpoints"""

        response = requests.get("/".join([self.explorer_api_url_base, "api/v0.3/s3_uri"]))
        response.raise_for_status()
        body = response.json()
        url_base = "/".join(
            [self.api_url_base, "cellxgene", "s3_uri", quote(quote(body, safe=""), safe=""), "api/v0.3"]
        )
        endpoints = ["config", "schema", "colors", "layout/obs", "annotations/obs", "annotations/var"]

        for endpoint in endpoints:
            with self.subTest(endpoint):
                # Test
                url = "/".join([url_base, endpoint])
                response = requests.head(url)

                # Verify
                response.raise_for_status()
                self.assertIn(response.headers["x-cache"], [self.cloudfront_miss, self.cloudfront_hit])

                # Call the request twice to make sure the cache wasn't cold.
                response = requests.head(url)
                response.raise_for_status()
                self.assertEqual(
                    response.headers["x-cache"],
                    self.cloudfront_hit,
                    msg=self.generate_error_msg(response, "Hit Expected"),
                )
