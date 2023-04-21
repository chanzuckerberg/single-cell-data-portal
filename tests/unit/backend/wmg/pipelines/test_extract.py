import unittest

import backend.wmg.pipeline.integrated_corpus.extract
from backend.wmg.data.constants import INCLUDED_ASSAYS


class TestExtract(unittest.TestCase):
    """
    Test case for extracting data to load into corpus to compute gene exression summary statistic cube
    """

    def generate_metadata(
        self, dataset_id, is_primary_data=True, include_h5ad=True, use_bad_assay=False, use_bad_organisms=False
    ):
        """
        Generate a set of dataset metadata that verifies the criteria for inclusion in the
        extraction pipeline
        """
        assay_ontologies = list(INCLUDED_ASSAYS.keys())
        metadata = dict()

        if is_primary_data:
            metadata["is_primary_data"] = [True]
        else:
            metadata["is_primary_data"] = [True, False]

        if use_bad_assay:
            metadata["assay"] = [dict(ontology_term_id=assay_ontologies[0], label="test_assay")]
        else:
            metadata["assay"] = [dict(ontology_term_id="bad_assay", label="test_bad_assay")]

        if use_bad_organisms:
            metadata["organism"] = [
                dict(ontology_term_id="NCBITaxon:9606", label="Homo Sapiens"),
                dict(ontology_term_id="NCBITaxon:Bad1", label="Bad1"),
                dict(ontology_term_id="NCBITaxon:Bad2", label="Bad2"),
            ]
        else:
            metadata["organism"] = [dict(ontology_term_id="NCBITaxon:9606", label="Homo Sapiens")]

        metadata["dataset_id"] = dataset_id

        dataset_assets = [
            {"filesize": 10, "filetype": "RDS", "url": "https://fake-url.rds"},
        ]
        if include_h5ad:
            dataset_assets.append({"filesize": 10, "filetype": "H5AD", "url": "https://fake-url.h5ad"})

        metadata["assets"] = dataset_assets
        return metadata

    def setUp(self):
        super().setUp()

        assay_ontologies = list(INCLUDED_ASSAYS.keys())

        # Dataset 1: should be included in the extraction pipeline
        self.dataset_0 = self.generate_metadata("dataset_0")

        # Dataset 2: should be included in the extraction pipeline
        self.dataset_1 = self.generate_metadata("dataset_1")
        # Only one assay needs to be included in the list of allowed assays
        self.dataset_1["assay"] = [
            dict(ontology_term_id="any_other_obo_id", label="test_assay"),
            dict(ontology_term_id=assay_ontologies[0], label="test_assay"),
        ]

        # Dataset 3: should NOT be included in the pipeline, since it has multiple organisms
        self.dataset__multiple_organisms = self.generate_metadata("dataset__multiple_organisms", use_bad_organisms=True)

        # Dataset 4: should NOT be included in the pipeline, since the assay isn't included in the whitelist
        self.dataset__wrong_assay = self.generate_metadata("dataset__wrong_assay", use_bad_assay=True)

        # Dataset 5: should NOT be included in the pipeline because is_primary_data is SECONDARY
        self.dataset__not_primary = self.generate_metadata("dataset__not_primary", is_primary_data=False)

        # Dataset 6: should NOT be included in the pipeline because it has no organism
        self.dataset__null_organism = self.generate_metadata("dataset__null_organism")
        self.dataset__null_organism["organism"] = None

    def test_get_asset_urls_pulls_expected_datasets(self):
        """
        Should only pull asset_urls for H5AD files
        Datasets should be
            - published
            - part of a public collection
            - contain primary data
            - not tombstoned
            - correct assay type
            - only contain one organism
        """
        expected_asset_urls = []
        not_expected_asset_urls = []
        expected_datasets = [self.dataset_0, self.dataset_1]
        not_expected_datasets = [
            self.dataset__not_primary,
            self.dataset__wrong_assay,
            self.dataset__multiple_organisms,
            self.dataset__null_organism,
        ]

        for dataset in not_expected_datasets:
            for asset in dataset["assets"]:
                not_expected_asset_urls.append(asset["url"])

        for dataset in expected_datasets:
            for asset in dataset["assets"]:
                if asset["filetype"] == "H5AD":
                    expected_asset_urls.append(asset["url"])

        asset_urls = backend.wmg.pipeline.integrated_corpus.extract.get_dataset_asset_urls(
            datasets=expected_datasets + not_expected_datasets
        ).values()

        self.assertCountEqual(asset_urls, expected_asset_urls)
        for uri in not_expected_datasets:
            self.assertNotIn(uri, asset_urls)
