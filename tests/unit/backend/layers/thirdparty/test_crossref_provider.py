import json
import unittest
from unittest.mock import patch

from requests import RequestException  # type: ignore
from requests.models import HTTPError, Response  # type: ignore

from backend.common.providers.crossref_provider import (
    CrossrefDOINotFoundException,
    CrossrefException,
    CrossrefFetchException,
    CrossrefParseException,
    CrossrefProvider,
)


class TestCrossrefProvider(unittest.TestCase):
    @patch("backend.common.providers.crossref_provider.requests.get")
    def test__provider_does_not_call_crossref_in_test(self, mock_get):
        provider = CrossrefProvider()
        res = provider.fetch_metadata("test_doi")
        self.assertIsNone(res)
        mock_get.assert_not_called()

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
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
                        "published-online": {"date-parts": [[2021, 11, 10]]},
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
            "published_at": 1636502400.0,
            "journal": "Nature",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, res)

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
    def test__provider_parses_authors_and_dates_correctly(self, mock_config, mock_get):

        response = Response()
        response.status_code = 200
        response._content = str.encode(
            json.dumps(
                {
                    "status": "ok",
                    "message": {
                        "author": [
                            {},
                            {"name": "A consortium"},
                            {"family": "Foo consortium"},
                            {
                                "family": "Smith",
                                "name": "Bar consortium",
                            },
                            {"given": "Jane"},
                            {
                                "given": "John",
                                "name": "Baz consortium",
                            },
                            {
                                "given": "John",
                                "family": "Doe",
                                "sequence": "first",
                            },
                            {
                                "given": "Jane",
                                "family": "Doe",
                                "name": "Bat consortium",
                            },
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
            "authors": [
                {"name": "A consortium"},
                {"name": "Foo consortium"},
                {"name": "Smith"},
                {"name": "Baz consortium"},
                {"given": "John", "family": "Doe"},
                {"given": "Jane", "family": "Doe"},
            ],
            "published_year": 2021,
            "published_month": 11,
            "published_day": 1,
            "published_at": 1635724800.0,
            "journal": "Nature",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, res)

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
    def test__provider_throws_exception_if_request_fails(self, mock_config, mock_get):
        """
        Asserts a CrossrefFetchException if the GET request fails for any reason
        """
        mock_get.side_effect = RequestException("Mocked CrossrefFetchException")

        provider = CrossrefProvider()

        with self.assertRaises(CrossrefFetchException):
            provider.fetch_metadata("test_doi")

        # Make sure that the parent CrossrefException will also be caught
        with self.assertRaises(CrossrefException):
            provider.fetch_metadata("test_doi")

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
    def test__provider_throws_exception_if_request_fails_with_404(self, mock_config, mock_get):
        """
        Asserts a CrossrefFetchException if the GET request fails for any reason
        """
        response_404 = Response()
        response_404.status_code = 404
        mock_get.side_effect = HTTPError(response=response_404)

        provider = CrossrefProvider()

        with self.assertRaises(CrossrefDOINotFoundException):
            provider.fetch_metadata("test_doi")

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
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

    @patch("backend.common.providers.crossref_provider.requests.get")
    @patch("backend.common.providers.crossref_provider.CorporaConfig")
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


if __name__ == "__main__":
    unittest.main()
