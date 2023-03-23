from dataclasses import asdict

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
        self.dataset_0 = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

        # Dataset 2: should be included in the extraction pipeline
        metadata = self.generate_base_metadata()
        # Only one assay needs to be included in the list of allowed assays
        metadata.assay = [
            OntologyTermId(ontology_term_id="any_other_obo_id", label="test_assay"),
            OntologyTermId(ontology_term_id=assay_ontologies[0], label="test_assay"),
        ]
        self.dataset_1 = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

        # Dataset 3: should NOT be included in the pipeline, since it has multiple organisms
        metadata = self.generate_base_metadata()
        metadata.organism = [
            OntologyTermId(ontology_term_id="NCBITaxon:9606", label="Homo Sapiens"),
            OntologyTermId(ontology_term_id="NCBITaxon:10090", label="Mus musculus"),
        ]
        self.dataset__multiple_organisms = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

        # Dataset 4: should NOT be included in the pipeline, since the assay isn't included in the whitelist
        metadata = self.generate_base_metadata()
        metadata.assay = [OntologyTermId(ontology_term_id="any_other_obo_id", label="test_assay")]
        self.dataset__wrong_assay = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

        # Dataset 5: should NOT be included in the pipeline because is_primary_data is SECONDARY
        metadata = self.generate_base_metadata()
        metadata.is_primary_data = "SECONDARY"
        self.dataset__not_primary = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

        # Dataset 6: should NOT be included in the pipeline because it doesn't belong to a published collection
        self.dataset__not_published = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(publish=False)
        )

        # Dataset 7: should NOT be included in the pipeline because it belongs to a tombstoned collection
        dataset_tmp = self.generate_dataset(publish=True)
        self.dataset__tombstoned = self.convert_business_logic_data_to_endpoint_dict_format(dataset_tmp)
        self.business_logic.tombstone_collection(CollectionId(dataset_tmp.collection_id))

        # Dataset 8: should NOT be included in the pipeline because it has no organism
        metadata = self.generate_base_metadata()
        metadata.organism = None
        self.dataset__null_organism = self.convert_business_logic_data_to_endpoint_dict_format(
            self.generate_dataset(metadata=metadata, publish=True)
        )

    def _convert_list_of_dataclass_to_list_of_dicts(self, list_of_dataclasses):
        """
        Convert a list of dataclasses to a list of dicts
        """
        if list_of_dataclasses is None:
            return None
        return [asdict(a) if a else a for a in list_of_dataclasses]

    def convert_business_logic_data_to_endpoint_dict_format(self, dataset):
        """
        Convert the data returned by the business logic into the format returned by the endpoint
        """
        dataset_business = self.business_logic.get_dataset_version(DatasetVersionId(dataset.dataset_version_id))
        dataset_assets = self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset.dataset_version_id))
        dataset_assets_list_of_dicts = []
        for asset in dataset_assets:
            asset = asdict(asset)
            asset_dict = {
                "filetype": asset["type"].upper(),
                "id": asset["id"]["id"],
                "s3_uri": asset["uri"],
                "filename": asset["uri"].split("/")[-1],
            }
            dataset_assets_list_of_dicts.append(asset_dict)

        dataset_dict = {}
        dataset_dict["dataset_version_id"] = dataset.dataset_version_id
        dataset_dict["dataset_id"] = dataset.dataset_id
        dataset_dict["collection_version_id"] = dataset.collection_version_id
        dataset_dict["collection_id"] = dataset.collection_id
        dataset_dict["artifact_ids"] = dataset.artifact_ids
        dataset_dict["is_primary_data"] = dataset_business.metadata.is_primary_data
        dataset_dict["organism"] = self._convert_list_of_dataclass_to_list_of_dicts(dataset_business.metadata.organism)
        dataset_dict["assay"] = self._convert_list_of_dataclass_to_list_of_dicts(dataset_business.metadata.assay)
        dataset_dict["explorer_url"] = dataset.explorer_url
        dataset_dict["dataset_assets"] = dataset_assets_list_of_dicts
        return dataset_dict

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
            dataset_assets = self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset["dataset_version_id"]))
            for asset in dataset_assets:
                not_expected_s3_uris.append(asset.uri)

        for dataset in expected_datasets:
            assets = self.business_logic.get_dataset_artifacts(DatasetVersionId(dataset["dataset_version_id"]))
            for asset in assets:
                if asset.type == DatasetArtifactType.H5AD:
                    expected_s3_uris.append(asset.uri)
        s3_uris = backend.wmg.pipeline.integrated_corpus.extract.get_dataset_s3_uris(
            datasets=expected_datasets + not_expected_datasets
        ).values()

        self.assertCountEqual(s3_uris, expected_s3_uris)
        for uri in not_expected_datasets:
            self.assertNotIn(uri, s3_uris)
