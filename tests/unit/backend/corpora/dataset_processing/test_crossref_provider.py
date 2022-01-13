import unittest
from unittest.mock import patch

from backend.corpora.common.crossref_provider import CrossrefProvider
import backend.corpora.common.crossref_provider


class TestCrossrefProvider(unittest.TestCase):
    @patch("backend.corpora.common.crossref_provider.requests.get")
    def test__provider_does_not_call_crossref_in_test(self, mock_get):
        provider = CrossrefProvider()
        res = provider.fetch_metadata("test_doi")
        self.assertIsNone(res)
        mock_get.assert_not_called()

    @patch("backend.corpora.common.crossref_provider.requests.get")
    def test__provider_calls_crossref_if_api_url_defined(self, mock_get):
        mock_get.return_value = {
            "status": "ok",
            "message": {
                "author": [
                    {
                        "given": "John",
                        "family": "Doe",
                        "sequence": "first",
                    },
                    {
                        "given": "Jane",
                        "family": "Doe",
                        "sequence": "additional",
                    },
                ],
                "published-online": {"date-parts": [[2021, 11, 18]]},
                "container-title": ["Nature"],
            },
        }
        # Pass a fake URL to the constructor. This will allow the provider to proceed and call the (mocked) Crossref API
        provider = CrossrefProvider("http://fake-api.uri/")
        res = provider.fetch_metadata("test_doi")
        mock_get.assert_called_once()

        expected_response = {
            "authors": [{"given": "John", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
            "published_year": 2021,
            "published_month": 11,
            "journal": "Nature",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, res)

    def test_split(self):
        self.assertTrue(False)
