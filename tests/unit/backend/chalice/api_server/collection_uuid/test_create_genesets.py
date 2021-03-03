import json
from backend.corpora.common.corpora_orm import CollectionVisibility, generate_uuid
from backend.corpora.common.entities.geneset import Geneset
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token


class TestGenesetCreation(BaseAuthAPITest):
    # Note for now the validation is done on the frontend so we are only testing with clean data
    # If we open up the api in the future we will need to validate and test for dirty data here as well

    def test_geneset_creation__ok(self):
        data = {
            "gene_sets": [
                [
                    "geneset1",
                    "this is geneset 1",
                    [
                        "gene1",
                        "gene 1 description",
                        {"randomProperty": "randomPropertyWords", "anotherProprerty": "moreWords"},
                    ],
                    ["gene2", "gene 2 description", {}],
                    ["gene3", "gene 3 description"],
                    ["gene4", "", {"randomProperty": "randomPropertyWords", "anotherProprerty": "moreWords"}],
                    ["gene5"],
                ],
                ["geneset2", "this is geneset2", ["onegene"]],
                [
                    "geneset3",
                    "this is geneset3",
                    ["1", "words", {"a": "b", "c": "d"}],
                    ["2", "words", {"a": "b", "c": "d"}],
                ],
            ]
        }
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        response.raise_for_status()

        body = json.loads(response.body)
        self.assertEqual(len(body), 3)
        geneset_names = [x["name"] for x in body]
        self.assertIn("geneset1", geneset_names)
        self.assertIn("geneset2", geneset_names)
        self.assertIn("geneset3", geneset_names)

    def test_geneset_creation__geneset_already_exists(self):
        data = {"gene_sets": [["geneset1", "this is geneset 1", ["gene1", "gene 1 description"]]]}
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        Geneset.create(self.session, collection=collection, name="geneset1", description="words", gene_symbols=["aaa"])
        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        self.assertEqual(400, response.status_code)

    def test_geneset_creation__not_owner_403(self):
        data = {"gene_sets": [["geneset1", "this is geneset 1", ["gene1", "gene 1 description"]]]}
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)

    def test_geneset_creation__not_private_403(self):
        data = {"gene_sets": [["geneset1", "this is geneset 1", ["gene1", "gene 1 description"]]]}
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)

    def test_geneset_creation__no_collection_403(self):
        uuid = generate_uuid()
        data = {"gene_sets": [["geneset1", "this is geneset 1", ["gene1", "gene 1 description"]]]}
        test_url = f"/dp/v1/collections/{uuid}/genesets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)
