import unittest
from unittest.mock import patch
from requests.models import Response
import json

from backend.corpora.common.providers.crossref_provider import (
    CrossrefException,
    CrossrefFetchException,
    CrossrefProvider,
    CrossrefParseException,
)


class TestCrossrefProvider(unittest.TestCase):
    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    def test__provider_does_not_call_crossref_in_test(self, mock_get):
        provider = CrossrefProvider()
        res = provider.fetch_metadata("test_doi")
        self.assertIsNone(res)
        mock_get.assert_not_called()

    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    @patch("backend.corpora.common.providers.crossref_provider.CorporaConfig")
    def test__provider_calls_crossref_if_api_key_defined(self, mock_config, mock_get):

        # Defining a mocked CorporaConfig will allow the provider to consider the `crossref_api_key`
        # not None, so it will go ahead and do the mocked call.

        response = Response()
        response.status_code = 200
        response._content = str.encode(
            json.dumps(
                {
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
            )
        )

        mock_get.return_value = response
        provider = CrossrefProvider()
        res = provider.fetch_metadata("test_doi")
        mock_get.assert_called_once()

        expected_response = {
            "authors": [{"given": "John", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
            "published_year": 2021,
            "published_month": 11,
            "published_day": 10,
            "published_at": 1636520400.0,
            "journal": "Nature",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, res)

    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    @patch("backend.corpora.common.providers.crossref_provider.CorporaConfig")
    def test__provider_parses_authors_and_dates_correctly(self, mock_config, mock_get):

        response = Response()
        response.status_code = 200
        response._content = str.encode(
            json.dumps(
                {
                    "status": "ok",
                    "message": {
                        "author": [
                            {
                                "given": "John",
                                "family": "Doe",
                                "sequence": "first",
                            },
                            {"name": "A consortium"},
                        ],
                        "published-online": {"date-parts": [[2021, 11]]},
                        "container-title": ["Nature"],
                    },
                }
            )
        )

        mock_get.return_value = response
        provider = CrossrefProvider()
        res = provider.fetch_metadata("test_doi")
        mock_get.assert_called_once()

        expected_response = {
            "authors": [{"given": "John", "family": "Doe"}, {"name": "A Consortium"}],
            "published_year": 2021,
            "published_month": 11,
            "published_day": 1,
            "published_at": 1635739200.0,
            "journal": "Nature",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, res)

    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    @patch("backend.corpora.common.providers.crossref_provider.CorporaConfig")
    def test__provider_throws_exception_if_request_fails(self, mock_config, mock_get):
        """
        Asserts a CrossrefFetchException if the GET request fails for any reason
        """
        mock_get.side_effect = Exception("Mocked CrossrefFetchException")

        provider = CrossrefProvider()

        with self.assertRaises(CrossrefFetchException):
            provider.fetch_metadata("test_doi")

        # Make sure that the parent CrossrefException will also be caught
        with self.assertRaises(CrossrefException):
            provider.fetch_metadata("test_doi")

    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    @patch("backend.corpora.common.providers.crossref_provider.CorporaConfig")
    def test__provider_throws_exception_if_request_fails_with_non_2xx_code(self, mock_config, mock_get):
        """
        Asserts a CrossrefFetchException if the GET request return a 500 error (any non 2xx will work)
        """

        response = Response()
        response.status_code = 500
        mock_get.return_value = response

        provider = CrossrefProvider()

        with self.assertRaises(CrossrefFetchException):
            provider.fetch_metadata("test_doi")

    @patch("backend.corpora.common.providers.crossref_provider.requests.get")
    @patch("backend.corpora.common.providers.crossref_provider.CorporaConfig")
    def test__provider_throws_exception_if_request_cannot_be_parsed(self, mock_config, mock_get):
        """
        Asserts an CrossrefParseException if the GET request succeeds but cannot be parsed
        """

        # Mocks a response that
        response = Response()
        response.status_code = 200
        response._content = str.encode(json.dumps({"status": "error"}))
        mock_get.return_value = response

        provider = CrossrefProvider()

        with self.assertRaises(CrossrefParseException):
            provider.fetch_metadata("test_doi")

        # Make sure that the parent CrossrefException will also be caught
        with self.assertRaises(CrossrefException):
            provider.fetch_metadata("test_doi")
