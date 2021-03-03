import json
from backend.corpora.common.corpora_orm import CollectionVisibility, generate_uuid
from backend.corpora.common.entities.geneset import Geneset
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token


class TestGenesetCreation(BaseAuthAPITest):
    # Note for now the validation is done on the frontend so we are only testing with clean data
    # If we open up the api in the future we will need to validate and test for dirty data here as well
    def setUp(self):
        super().setUp()
        self.geneset1 = {
            "gene_set_name": "geneset1",
            "gene_set_description": "this is geneset 1",
            "genes": [
                {
                    "gene_symbol": "gene1",
                    "gene_description": "describe a gene",
                    "additional_params": {"randomProperty": "randomPropertyWords", "anotherProprerty": "moreWords"},
                },
                {"gene_symbol": "1", "gene_description": "words", "additional_params": {"a": "b", "c": "d"}},
            ],
        }
        self.geneset2 = {
            "gene_set_name": "geneset2",
            "gene_set_description": "this is geneset 2",
            "genes": [{"gene_symbol": "onegene"}],
        }
        self.geneset3 = {
            "gene_set_name": "geneset3",
            "gene_set_description": "this is geneset 3",
            "genes": [
                {"gene_symbol": "1", "gene_description": "words", "additional_params": {"a": "b", "c": "d"}},
                {"gene_symbol": "2", "gene_description": "words", "additional_params": {"a": "b", "c": "d"}},
                {"gene_symbol": "3", "gene_description": "words", "additional_params": {"a": "b", "c": "d"}},
                {"gene_symbol": "4", "gene_description": "words", "additional_params": {"a": "b", "c": "d"}},
            ],
        }
        self.collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        self.test_url = f"/dp/v1/collections/{self.collection.id}/genesets"

    def test_geneset_creation__ok(self):
        data = {"gene_sets": [self.geneset1, self.geneset2, self.geneset3]}
        response = self.app.post(self.test_url, headers=self.headers, data=json.dumps(data))
        response.raise_for_status()

        body = json.loads(response.body)
        self.assertEqual(len(body), 3)
        geneset_names = [x["name"] for x in body]
        self.assertIn("geneset1", geneset_names)
        self.assertIn("geneset2", geneset_names)
        self.assertIn("geneset3", geneset_names)

    def test_all_gene_information_correctly_stored(self):
        data = {"gene_sets": [self.geneset1, self.geneset2, self.geneset3]}
        response = self.app.post(self.test_url, headers=self.headers, data=json.dumps(data))
        response.raise_for_status()
        body = json.loads(response.body)
        genesets = {}
        for x in body:
            genesets[x["name"]] = x["id"]
        geneset1 = Geneset.get(self.session, genesets["geneset1"])
        geneset2 = Geneset.get(self.session, genesets["geneset2"])
        geneset3 = Geneset.get(self.session, genesets["geneset3"])

        self.assertEqual(geneset1.gene_symbols, self.geneset1["genes"])
        self.assertEqual(geneset2.gene_symbols, self.geneset2["genes"])
        self.assertEqual(geneset3.gene_symbols, self.geneset3["genes"])

    def test_geneset_creation__geneset_already_exists(self):
        data = {"gene_sets": [self.geneset1]}
        Geneset.create(
            self.session, collection=self.collection, name="geneset1", description="words", gene_symbols=["aaa"]
        )
        response = self.app.post(self.test_url, headers=self.headers, data=json.dumps(data))
        self.assertEqual(400, response.status_code)

    def test_geneset_creation__not_owner_403(self):
        data = {"gene_sets": [self.geneset1, self.geneset2]}
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"

        response = self.app.post(test_url, headers=self.headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)

    def test_geneset_creation__not_private_403(self):
        data = {"gene_sets": [self.geneset2, self.geneset3]}
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        test_url = f"/dp/v1/collections/{collection.id}/genesets"

        response = self.app.post(test_url, headers=self.headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)

    def test_geneset_creation__no_collection_403(self):
        uuid = generate_uuid()
        data = {"gene_sets": [self.geneset1]}
        test_url = f"/dp/v1/collections/{uuid}/genesets"
        response = self.app.post(test_url, headers=self.headers, data=json.dumps(data))
        self.assertEqual(403, response.status_code)
