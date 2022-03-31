import unittest

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.wmg.data import extract
from backend.wmg.data.extract import included_assay_ontologies
from tests.unit.backend.fixtures.generate_data_mixin import GenerateDataMixin
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestExtract(CorporaTestCaseUsingMockAWS, GenerateDataMixin):
    """
    Test case for extracting data to load into corpus to compute gene exression summary statistic cube
    """

    def setUp(self):
        super().setUp()
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        # INCLUDE
        assay_ontologies = list(included_assay_ontologies.keys())
        self.dataset_0 = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tombstone=False,
            assay={"ontology_term_id": assay_ontologies[0], "label": "test_assay"},
        )
        self.dataset_1 = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tombstone=False,
            assay={"ontology_term_id": assay_ontologies[1], "label": "test_assay"},
        )

        # DONT INCLUDE
        self.dataset__wrong_assay = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tombstone=False,
            assay={"ontology_term_id": "any_other_obo_id", "label": "test_assay"},
        )
        self.dataset__not_primary = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="SECONDARY",
            tombstone=True,
            assay={"ontology_term_id": assay_ontologies[1], "label": "test_assay"},
        )

        self.dataset__not_published = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=False,
            is_primary_data="PRIMARY",
            tombstone=False,
            assay={"ontology_term_id": assay_ontologies[1], "label": "test_assay"},
        )
        self.dataset__tombstoned = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PUBLIC",
            collection_id=pub_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tombstone=True,
            assay={"ontology_term_id": assay_ontologies[1], "label": "test_assay"},
        )

        private_collection = self.generate_collection(self.session, visibility="PRIVATE")
        self.dataset__private_collection = self.generate_dataset_with_s3_resources(
            self.session,
            artifacts=True,
            explorer_s3_object=False,
            collection_visibility="PRIVATE",
            collection_id=private_collection.id,
            published=True,
            is_primary_data="PRIMARY",
            tombstone=False,
            assay={"ontology_term_id": assay_ontologies[0], "label": "test_assay"},
        )

    def tearDown(self):
        super().tearDown()

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
        expected_s3_uris = []
        not_expected_s3_uris = []
        expected_datasets = [self.dataset_0, self.dataset_1]
        not_expected_datasets = [
            self.dataset__tombstoned,
            self.dataset__not_primary,
            self.dataset__not_published,
            self.dataset__private_collection,
            self.dataset__wrong_assay,
        ]
        for dataset in not_expected_datasets:
            dataset_assets = dataset.get_assets()
            for asset in dataset_assets:
                not_expected_s3_uris.append(asset.s3_uri)
        for dataset in expected_datasets:
            assets = dataset.get_assets()
            for asset in assets:
                if asset.filetype == DatasetArtifactFileType.H5AD:
                    expected_s3_uris.append(asset.s3_uri)

        s3_uris = set(extract.get_dataset_s3_uris())
        self.assertEquals(set(expected_s3_uris), s3_uris)

        @unittest.skip
        def test_datasets_copied_to_correct_location(self):
            # TODO implement with changes from @ebezzis pipeline work
            raise NotImplementedError
