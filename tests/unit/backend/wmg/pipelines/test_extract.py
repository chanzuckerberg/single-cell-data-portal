import backend.wmg.pipeline.integrated_corpus.extract
from backend.layers.common.entities import CollectionId, DatasetArtifactType, DatasetVersionId, OntologyTermId
from backend.wmg.data.constants import INCLUDED_ASSAYS
from tests.unit.backend.layers.common.base_test import BaseTest


class TestExtract(BaseTest):
    """
    Test case for extracting data to load into corpus to compute gene exression summary statistic cube
    """

    def generate_base_metadata(self):
        """
        Generate a set of dataset metadata that verifies the criteria for inclusion in the
        extraction pipeline
        """
        assay_ontologies = list(INCLUDED_ASSAYS.keys())
        metadata = self.get_sample_dataset_metadata()
        metadata.is_primary_data = "PRIMARY"
        metadata.assay = [OntologyTermId(ontology_term_id=assay_ontologies[0], label="test_assay")]
        metadata.organism = [OntologyTermId(ontology_term_id="NCBITaxon:9606", label="Homo Sapiens")]
        return metadata

    def setUp(self):
        super().setUp()

        assay_ontologies = list(INCLUDED_ASSAYS.keys())

        # Dataset 1: should be included in the extraction pipeline
        metadata = self.generate_base_metadata()
        self.dataset_0 = self.generate_dataset(metadata=metadata, publish=True)

        # Dataset 2: should be included in the extraction pipeline
        metadata = self.generate_base_metadata()
        # Only one assay needs to be included in the list of allowed assays
        metadata.assay = [
            OntologyTermId(ontology_term_id="any_other_obo_id", label="test_assay"),
            OntologyTermId(ontology_term_id=assay_ontologies[0], label="test_assay"),
        ]
        self.dataset_1 = self.generate_dataset(metadata=metadata, publish=True)

        # Dataset 3: should NOT be included in the pipeline, since it has multiple organisms
        metadata = self.generate_base_metadata()
        metadata.organism = [
            OntologyTermId(ontology_term_id="NCBITaxon:9606", label="Homo Sapiens"),
            OntologyTermId(ontology_term_id="NCBITaxon:10090", label="Mus musculus"),
        ]
        self.dataset__multiple_organisms = self.generate_dataset(metadata=metadata, publish=True)

        # Dataset 4: should NOT be included in the pipeline, since the assay isn't included in the whitelist
        metadata = self.generate_base_metadata()
        metadata.assay = [OntologyTermId(ontology_term_id="any_other_obo_id", label="test_assay")]
        self.dataset__wrong_assay = self.generate_dataset(metadata=metadata, publish=True)

        # Dataset 5: should NOT be included in the pipeline because is_primary_data is SECONDARY
        metadata = self.generate_base_metadata()
        metadata.is_primary_data = "SECONDARY"
        self.dataset__not_primary = self.generate_dataset(metadata=metadata, publish=True)

        # Dataset 6: should NOT be included in the pipeline because it doesn't belong to a published collection
        self.dataset__not_published = self.generate_dataset(publish=False)

        # Dataset 7: should NOT be included in the pipeline because it belongs to a tombstoned collection
        self.dataset__tombstoned = self.generate_dataset(publish=True)
        self.business_logic.tombstone_collection(CollectionId(self.dataset__tombstoned.collection_id))

        # Dataset 8: should NOT be included in the pipeline because it has no organism
        metadata = self.generate_base_metadata()
        metadata.organism = None
        self.dataset__null_organism = self.generate_dataset(metadata=metadata, publish=True)

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
            - only contain one organism
        """
        expected_s3_uris = []
        not_expected_s3_uris = []
        expected_datasets = [self.dataset_0, self.dataset_1]
        not_expected_datasets = [
            self.dataset__tombstoned,
            self.dataset__not_primary,
            self.dataset__not_published,
            self.dataset__wrong_assay,
            self.dataset__multiple_organisms,
            self.dataset__null_organism,
        ]

        for dataset in not_expected_datasets:
            dataset_assets = self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset.dataset_version_id))
            for asset in dataset_assets:
                not_expected_s3_uris.append(asset.uri)

        for dataset in expected_datasets:
            assets = self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset.dataset_version_id))
            for asset in assets:
                if asset.type == DatasetArtifactType.H5AD:
                    expected_s3_uris.append(asset.uri)

        s3_uris = backend.wmg.pipeline.integrated_corpus.extract.get_dataset_s3_uris().values()

        self.assertCountEqual(s3_uris, expected_s3_uris)
        for uri in not_expected_datasets:
            self.assertNotIn(uri, s3_uris)
