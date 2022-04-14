import json

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities.geneset import GenesetDatasetLink
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest, BasicAuthAPITestCurator
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestGenesets(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def _get_geneset_ids(self, collection_id, headers, public=False):
        if public:
            rsp = self.app.get(f"/dp/v1/collections/{collection_id}", headers=headers)
        else:
            rsp = self.app.get(f"/dp/v1/collections/{collection_id}?visibility=PRIVATE", headers=headers)
        self.assertEqual(200, rsp.status_code)
        bdy = json.loads(rsp.data)
        return [g["id"] for g in bdy.get("genesets", [])]

    def _delete_geneset_test(self, collection_id, headers, geneset):
        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection_id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)

        # delete the geneset
        actual_geneset_ids = self._get_geneset_ids(collection_id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)
        response = self.app.delete(f"/dp/v1/genesets/{geneset.id}", headers=headers)
        self.assertEqual(202, response.status_code)

        # try to get geneset
        actual_geneset_ids = self._get_geneset_ids(collection_id, headers)
        self.assertNotIn(geneset.id, actual_geneset_ids)

        # Delete a seconds time
        response = self.app.delete(f"/dp/v1/genesets/{geneset.id}", headers=headers)
        self.assertEqual(202, response.status_code)

    def test__delete_geneset_ACCEPTED(self):
        """
        Delete a geneset from a private collection.
        """

        headers = dict(host="localhost", Cookie=get_auth_token(self.app))
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE, owner="test_user_id"
        )

        with self.subTest("With collection"):
            # create a geneset
            geneset = self.generate_geneset(
                self.session, collection_id=collection.id, collection_visibility=collection.visibility
            )

            # get geneset
            actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
            self.assertIn(geneset.id, actual_geneset_ids)

            self._delete_geneset_test(collection.id, headers, geneset)

        with self.subTest("With dataset"):
            # create dataset
            dataset = self.generate_dataset(
                self.session, collection_id=collection.id, collection_visibility=collection.visibility
            )

            # create a geneset
            geneset = self.generate_geneset(
                self.session, collection_id=collection.id, collection_visibility=collection.visibility
            )
            GenesetDatasetLink.create(self.session, geneset.id, dataset.id)

            # get geneset
            actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
            self.assertIn(geneset.id, actual_geneset_ids)

            self._delete_geneset_test(collection.id, headers, geneset)

            # get dataset
            response = self.app.get(f"/dp/v1/collections/{collection.id}?visibility=PRIVATE", headers=headers)
            self.assertEqual(200, response.status_code)
            body = json.loads(response.data)
            self.assertIn(dataset.id, [d["id"] for d in body.get("datasets", [])])

    def test__delete_gene_set__UNAUTHORIZED(self):
        headers = dict(host="localhost", Cookie=get_auth_token(self.app))
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE, owner="some_one_else"
        )

        # create a geneset
        geneset = self.generate_geneset(
            self.session, collection_id=collection.id, collection_visibility=collection.visibility
        )

        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)

        # delete the geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)
        response = self.app.delete(f"/dp/v1/genesets/{geneset.id}", headers=headers)
        self.assertEqual(403, response.status_code)

        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)

    def test__delete_gene_set__NOT_ALLOWED(self):
        headers = dict(host="localhost", Cookie=get_auth_token(self.app))
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC, owner="test_user_id"
        )

        # create a geneset
        geneset = self.generate_geneset(
            self.session, collection_id=collection.id, collection_visibility=collection.visibility
        )

        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers, True)
        self.assertIn(geneset.id, actual_geneset_ids)

        # delete the geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers, True)
        self.assertIn(geneset.id, actual_geneset_ids)
        response = self.app.delete(f"/dp/v1/genesets/{geneset.id}", headers=headers)
        self.assertEqual(405, response.status_code)

        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers, True)
        self.assertIn(geneset.id, actual_geneset_ids)


class TestGenesetsCurators(BasicAuthAPITestCurator, CorporaTestCaseUsingMockAWS):
    def _get_geneset_ids(self, collection_id, headers, public=False):
        if public:
            rsp = self.app.get(f"/dp/v1/collections/{collection_id}", headers=headers)
        else:
            rsp = self.app.get(f"/dp/v1/collections/{collection_id}?visibility=PRIVATE", headers=headers)
        self.assertEqual(200, rsp.status_code)
        bdy = json.loads(rsp.data)
        return [g["id"] for g in bdy.get("genesets", [])]

    def test__delete_gene_set__200(self):
        headers = dict(host="localhost", Cookie=get_auth_token(self.app))
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE, owner="some_one_else"
        )

        # create a geneset
        geneset = self.generate_geneset(
            self.session, collection_id=collection.id, collection_visibility=collection.visibility
        )

        # get geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)

        # delete the geneset
        actual_geneset_ids = self._get_geneset_ids(collection.id, headers)
        self.assertIn(geneset.id, actual_geneset_ids)
        response = self.app.delete(f"/dp/v1/genesets/{geneset.id}", headers=headers)
        self.assertEqual(202, response.status_code)
