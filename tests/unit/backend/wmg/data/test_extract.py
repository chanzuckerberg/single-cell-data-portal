import unittest

from backend.corpora.common.entities import DatasetAsset
from backend.wmg.data import extract
from backend.wmg.data.extract import included_assay_ontology_ids
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestExtract(CorporaTestCaseUsingMockAWS):
    """
    Test case for extracting data to load into corpus to compute gene exression summary statistic cube
    """

    def setUp(self):
        super().setUp()
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        ## INCLUDE
        self.dataset_0 = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tomstoned="f",
            assay={"ontology_term_id": included_assay_ontology_ids[0], "label": "test_assay"}
        )
        self.dataset_1 = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tomstoned="f",
            assay={"ontology_term_id": included_assay_ontology_ids[1], "label": "test_assay"}
        )

        ## DONT INCLUDE
        self.dataset__wrong_assay = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tomstoned="f",
            assay={"ontology_term_id": "any_other_obo_id", "label": "test_assay"}
        )
        self.dataset__not_primary_0 = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="SECONDARY",
            tomstoned="t",
            assay={"ontology_term_id": included_assay_ontology_ids[1], "label": "test_assay"}
        )
        self.dataset__not_primary_1 = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="both",
            tomstoned="t",
            assay={"ontology_term_id": included_assay_ontology_ids[1], "label": "test_assay"}
        )
        self.dataset__not_published = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=False,
            is_primary_data="PRIMARY",
            tomstoned="t",
            assay={"ontology_term_id": included_assay_ontology_ids[1], "label": "test_assay"}
        )
        self.dataset__tombstoned = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tomstoned="t",
            assay={"ontology_term_id": included_assay_ontology_ids[1], "label": "test_assay"}
        )

        private_collection = self.generate_collection(self.session, visibility="PRIVATE")
        self.dataset__private_collection = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PUBLIC",
            collection_id=private_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tomstoned="f",
            assay={"ontology_term_id": included_assay_ontology_ids[0], "label": "test_assay"}
        )

    def test_get_s3_uris_pulls_expected_datasets(self):
        """
        Should only pull s3_uris for H5AD files
        Datasets should be
            - published
            - part of a public collection
            - contain primary data
            - not tombstoned
            - correct assay type
        """
        s3_uris = extract.get_s3_uris()
        s3_uri_0 = self.dataset_0.get_assets()
        import pdb
        pdb.set_trace()
        self.assertEqual(1, 2)



        @unittest.skip
        def test_datasets_copied_to_correct_location(self):
            # TODO implement when local stack working correctly
            raise NotImplementedError
