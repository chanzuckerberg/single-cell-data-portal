import copy
import json
import unittest
from unittest.mock import Mock, patch

from requests import RequestException
from requests.models import HTTPError, Response

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
        metadata, doi = provider.fetch_metadata("test_doi")
        self.assertIsNone(metadata)
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
        res, doi_curie_from_crossref = provider.fetch_metadata("test_doi")
        self.assertEqual("test_doi", doi_curie_from_crossref)
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
    def test__published_doi_used_if_exists_for_preprint(self, mock_config, mock_get):

        # Defining a mocked CorporaConfig will allow the provider to consider the `crossref_api_key`
        # not None, so it will go ahead and do the mocked call.

        def make_response(content):
            response = Response()
            response.status_code = 200
            response._content = str.encode(json.dumps(content))
            return response

        body = {
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
                "subtype": "preprint",
                "relation": {
                    "is-preprint-of": [
                        {"id-type": "not_doi"},
                        {
                            "id": "published_doi",
                            "id-type": "doi",
                        },
                    ]
                },
            },
        }

        provider = CrossrefProvider()

        with self.subTest("Published DOI is used when available"):
            preprint_body = copy.deepcopy(body)
            response_preprint = make_response(preprint_body)

            published_body = copy.deepcopy(body)
            del published_body["message"]["subtype"]
            published_body["message"]["author"][0]["given"] = "Jonathan"
            response_published = make_response(published_body)

            responses = [response_published, response_preprint]
            mock_get.side_effect = lambda *x, **y: responses.pop()
            res, doi_curie_from_crossref = provider.fetch_metadata("preprint_doi")
            self.assertEqual("published_doi", doi_curie_from_crossref)
            expected_response = {
                "authors": [{"given": "Jonathan", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
                "published_year": 2021,
                "published_month": 11,
                "published_day": 10,
                "published_at": 1636502400.0,
                "journal": "Nature",
                "is_preprint": False,
            }

            self.assertDictEqual(expected_response, res)

        with (
            self.subTest("Preprint DOI is used when published is referenced but cannot be retrieved")
            and patch.object(provider, "fetch_published_metadata") as fetch_published_metadata_mock
        ):
            fetch_published_metadata_mock.return_value = (None, None)

            preprint_body = copy.deepcopy(body)
            response_preprint = make_response(preprint_body)

            published_body = copy.deepcopy(body)
            del published_body["message"]["subtype"]
            response_published = make_response(published_body)

            responses = [response_published, response_preprint]
            mock_get.side_effect = lambda *x, **y: responses.pop()
            res, doi_curie_from_crossref = provider.fetch_metadata("preprint_doi")
            self.assertEqual("preprint_doi", doi_curie_from_crossref)
            expected_response = {
                "authors": [{"given": "John", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
                "published_year": 2021,
                "published_month": 11,
                "published_day": 10,
                "published_at": 1636502400.0,
                "journal": "Nature",
                "is_preprint": True,
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
                            # Test case for dupe removal
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
        res, _ = provider.fetch_metadata("test_doi")
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
    def test__provider_unescapes_journal_correctly(self, mock_config, mock_get):
        response = Response()
        response.status_code = 200
        response._content = str.encode(
            json.dumps(
                {
                    "status": "ok",
                    "message": {
                        "author": [
                            {"name": "A consortium"},
                        ],
                        "published-online": {"date-parts": [[2021, 11]]},
                        "container-title": ["Clinical &amp; Translational Med"],
                    },
                }
            )
        )

        mock_get.return_value = response
        provider = CrossrefProvider()
        author_data, _ = provider.fetch_metadata("test_doi")
        mock_get.assert_called_once()

        expected_response = {
            "authors": [
                {"name": "A consortium"},
            ],
            "published_year": 2021,
            "published_month": 11,
            "published_day": 1,
            "published_at": 1635724800.0,
            "journal": "Clinical & Translational Med",
            "is_preprint": False,
        }

        self.assertDictEqual(expected_response, author_data)

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

    @patch("backend.common.providers.crossref_provider.requests.get")
    def test__get_title_and_citation_from_doi(self, mock_get):
        provider = CrossrefProvider()
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "message": {
                "title": ["Test Title"],
                "author": [{"family": "Doe", "given": "John"}],
                "container-title": ["Test Journal"],
                "created": {"date-parts": [[2022]]},
            }
        }
        mock_get.return_value = mock_response

        result = provider.get_title_and_citation_from_doi("10.1016/j.cell.2019.11.025")
        self.assertEqual(result, "Test Title\n\n - Doe et al. (2022) Test Journal")


if __name__ == "__main__":
    unittest.main()
