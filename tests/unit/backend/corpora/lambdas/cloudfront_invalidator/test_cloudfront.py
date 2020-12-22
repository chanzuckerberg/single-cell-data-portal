from unittest import mock, TestCase

from backend.corpora.lambdas.cloudfront_invalidator.cloudfront import (
    get_cloudfront_distribution,
    invalidate_distributions,
)


class TestCloudfront(TestCase):
    @mock.patch("backend.corpora.lambdas.cloudfront_invalidator.cloudfront._cloudfront_client")
    def test_get_cloudfront_distribution(self, mock_cloudfront_client):
        mock_cloudfront_client.list_distributions.return_value = {
            "DistributionList": {
                "Items": [
                    {"Id": "id1", "Origins": {"Quantity": 1, "Items": [{"DomainName": "test_bucket-one-match"}]}},
                    {"Id": "id2", "Origins": {"Quantity": 1, "Items": [{"DomainName": "test_bucket"}]}},
                ]
            }
        }
        with self.subTest("Single distribution"):
            self.assertListEqual(get_cloudfront_distribution("test_bucket-one-match"), ["id1"])
        with self.subTest("Multiple distributions"):
            self.assertListEqual(get_cloudfront_distribution("test_bucket"), ["id1", "id2"])
        with self.subTest("No Match"):
            self.assertListEqual(get_cloudfront_distribution("no-match"), [])

    @mock.patch("backend.corpora.lambdas.cloudfront_invalidator.cloudfront._cloudfront_client")
    def test_invalid_distributions(self, mock_cloudfront_client):
        with self.subTest("No distributions"):
            mock_cloudfront_client.invalidate_distribution.return_value = []
            self.assertListEqual(invalidate_distributions([]), [])
        with self.subTest("Multiple Distributions"):
            distributions = ["distribution_1", "distribution_2"]
            mock_cloudfront_client.invalidate_distribution.return_value = [
                {"Invalidation": {"Id": d}} for d in distributions
            ]
            self.assertEqual(len(invalidate_distributions(distributions)), 2)
