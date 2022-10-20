import time
from unittest.mock import patch

import anndata

import numpy as np
import pandas
from moto import mock_s3

from backend.common.corpora_orm import UploadStatus
from backend.common.entities.dataset import Dataset
from backend.dataset_pipeline.processing.exceptions import ProcessingCancelled
from backend.dataset_pipeline.processing.process_download_validate import download_from_source_uri, extract_metadata
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestProcessingDownloadValidate(DataPortalTestCase):
    @patch("scanpy.read_h5ad")
    def test_extract_metadata(self, mock_read_h5ad):

        df = pandas.DataFrame(
            np.random.randint(10, size=(50001, 5)) * 50, columns=list("ABCDE"), index=(str(i) for i in range(50001))
        )

        self_reported_ethnicity = tissue = np.random.choice([0, 1], size=(50001))
        assay = development_stage = sex = np.random.choice([0, 1, 2], size=(50001))

        obs = pandas.DataFrame(
            np.hstack(
                [
                    np.array([["lung", "liver"][i] for i in tissue]).reshape(50001, 1),
                    np.array([["UBERON:01", "UBERON:10"][i] for i in tissue]).reshape(50001, 1),
                    np.array([["10x", "smartseq", "cite-seq"][i] for i in assay]).reshape(50001, 1),
                    np.array([["EFO:001", "EFO:010", "EFO:011"][i] for i in assay]).reshape(50001, 1),
                    np.random.choice(["healthy"], size=(50001, 1)),
                    np.random.choice(["MONDO:123"], size=(50001, 1)),
                    np.array([["male", "female", "fixed"][i] for i in sex]).reshape(50001, 1),
                    np.array([["M", "F", "MF"][i] for i in sex]).reshape(50001, 1),
                    np.array([["solomon islander", "orcadian"][i] for i in self_reported_ethnicity]).reshape(50001, 1),
                    np.array([["HANCESTRO:321", "HANCESTRO:456"][i] for i in self_reported_ethnicity]).reshape(
                        50001, 1
                    ),
                    np.array([["adult", "baby", "tween"][i] for i in development_stage]).reshape(50001, 1),
                    np.array([["HsapDv:0", "HsapDv:1", "HsapDv:2"][i] for i in development_stage]).reshape(50001, 1),
                    np.random.choice(["Homo sapiens"], size=(50001, 1)),
                    np.random.choice(["NCBITaxon:8505"], size=(50001, 1)),
                    np.random.choice([0], size=(50001, 1)),
                    np.random.choice(["liver"], size=(50001, 1)),
                    np.random.choice(["Hepatic-1A"], size=(50001, 1)),
                    np.array([["cell", "nucleus", "na"][i] for i in assay]).reshape(50001, 1),
                    np.random.choice(["F1", "F2"], size=(50001, 1)),
                ]
            ),
            columns=[
                "tissue",
                "tissue_ontology_term_id",
                "assay",
                "assay_ontology_term_id",
                "disease",
                "disease_ontology_term_id",
                "sex",
                "sex_ontology_term_id",
                "self_reported_ethnicity",
                "self_reported_ethnicity_ontology_term_id",
                "development_stage",
                "development_stage_ontology_term_id",
                "organism",
                "organism_ontology_term_id",
                "is_primary_data",
                "cell_type",
                "cell_type_ontology_term_id",
                "suspension_type",
                "donor_id",
            ],
            index=(str(i) for i in range(50001)),
        )
        uns = {
            "title": "my test dataset",
            "X_approximate_distribution": "normal",
            "batch_condition": np.array({"batchA", "batchB"}),
            "schema_version": "3.0.0",
        }

        var = pandas.DataFrame(
            data=["gene", "spike-in", "gene", "gene", "gene"], columns=["feature_biotype"], index=df.columns
        )

        adata = anndata.AnnData(X=df, obs=obs, uns=uns, var=var)
        mock_read_h5ad.return_value = adata

        extracted_metadata = extract_metadata("dummy")
        lab, ont = "label", "ontology_term_id"

        def list_equal(list1, list2, cmp_func):
            self.assertEqual(len(list1), len(list2))
            for el1 in list1:
                self.assertIn(el1, list2)
                el2 = list2[list2.index(el1)]
                if cmp_func:
                    cmp_func(el1, el2)

        list_equal(extracted_metadata["organism"], [{lab: "Homo sapiens", ont: "NCBITaxon:8505"}], self.assertDictEqual)

        list_equal(
            extracted_metadata["tissue"],
            [{lab: "lung", ont: "UBERON:01"}, {lab: "liver", ont: "UBERON:10"}],
            self.assertDictEqual,
        )

        list_equal(
            extracted_metadata["assay"],
            [{lab: "10x", ont: "EFO:001"}, {lab: "smartseq", ont: "EFO:010"}, {lab: "cite-seq", ont: "EFO:011"}],
            self.assertDictEqual,
        )

        list_equal(extracted_metadata["disease"], [{lab: "healthy", ont: "MONDO:123"}], self.assertDictEqual)

        list_equal(
            extracted_metadata["sex"],
            [{lab: "male", ont: "M"}, {lab: "female", ont: "F"}, {lab: "fixed", ont: "MF"}],
            self.assertDictEqual,
        )

        list_equal(
            extracted_metadata["self_reported_ethnicity"],
            [{lab: "solomon islander", ont: "HANCESTRO:321"}, {lab: "orcadian", ont: "HANCESTRO:456"}],
            self.assertDictEqual,
        )

        list_equal(
            extracted_metadata["development_stage"],
            [{lab: "adult", ont: "HsapDv:0"}, {lab: "baby", ont: "HsapDv:1"}, {lab: "tween", ont: "HsapDv:2"}],
            self.assertDictEqual,
        )

        list_equal(
            extracted_metadata["suspension_type"],
            ["cell", "nucleus", "na"],
            self.assertEqual,
        )

        list_equal(
            extracted_metadata["donor_id"],
            ["F1", "F2"],
            self.assertEqual,
        )

        self.assertEqual(extracted_metadata["x_approximate_distribution"], "NORMAL")
        self.assertEqual(extracted_metadata["batch_condition"], np.array({"batchA", "batchB"}))
        self.assertEqual(extracted_metadata["schema_version"], "3.0.0")

        self.assertEqual(extracted_metadata["cell_count"], 50001)

        filter = np.where(adata.var.feature_biotype == "gene")[0]
        self.assertAlmostEqual(extracted_metadata["mean_genes_per_cell"], np.count_nonzero(adata.X[:, filter]) / 50001)

    @patch("scanpy.read_h5ad")
    def test_extract_metadata_find_raw_layer(self, mock_read_h5ad):
        # Setup anndata to be read
        non_zeros_X_layer_df = pandas.DataFrame(
            np.full((11, 3), 2), columns=list("ABC"), index=(str(i) for i in range(11))
        )
        zeros_layer_df = pandas.DataFrame(np.zeros((11, 3)), columns=list("ABC"), index=(str(i) for i in range(11)))

        self_reported_ethnicity = tissue = np.random.choice([0, 1], size=(11))
        assay = development_stage = sex = np.random.choice([0, 1, 2], size=(11))

        obs = pandas.DataFrame(
            np.hstack(
                [
                    np.array([["lung", "liver"][i] for i in tissue]).reshape(11, 1),
                    np.array([["UBERON:01", "UBERON:10"][i] for i in tissue]).reshape(11, 1),
                    np.array([["10x", "smartseq", "cite-seq"][i] for i in assay]).reshape(11, 1),
                    np.array([["EFO:001", "EFO:010", "EFO:011"][i] for i in assay]).reshape(11, 1),
                    np.random.choice(["healthy"], size=(11, 1)),
                    np.random.choice(["MONDO:123"], size=(11, 1)),
                    np.array([["male", "female", "fixed"][i] for i in sex]).reshape(11, 1),
                    np.array([["M", "F", "MF"][i] for i in sex]).reshape(11, 1),
                    np.array([["solomon islander", "orcadian"][i] for i in self_reported_ethnicity]).reshape(11, 1),
                    np.array([["HANCESTRO:321", "HANCESTRO:456"][i] for i in self_reported_ethnicity]).reshape(11, 1),
                    np.array([["adult", "baby", "tween"][i] for i in development_stage]).reshape(11, 1),
                    np.array([["HsapDv:0", "HsapDv:1", "HsapDv:2"][i] for i in development_stage]).reshape(11, 1),
                    np.random.choice(["Homo sapiens"], size=(11, 1)),
                    np.random.choice(["NCBITaxon:8505"], size=(11, 1)),
                    np.random.choice([0], size=(11, 1)),
                    np.random.choice(["liver"], size=(11, 1)),
                    np.random.choice(["Hepatic-1A"], size=(11, 1)),
                    np.array([["cell", "nucleus", "na"][i] for i in assay]).reshape(11, 1),
                    np.random.choice(["F1", "F2"], size=(11, 1)),
                ]
            ),
            columns=[
                "tissue",
                "tissue_ontology_term_id",
                "assay",
                "assay_ontology_term_id",
                "disease",
                "disease_ontology_term_id",
                "sex",
                "sex_ontology_term_id",
                "self_reported_ethnicity",
                "self_reported_ethnicity_ontology_term_id",
                "development_stage",
                "development_stage_ontology_term_id",
                "organism",
                "organism_ontology_term_id",
                "is_primary_data",
                "cell_type",
                "cell_type_ontology_term_id",
                "suspension_type",
                "donor_id",
            ],
            index=(str(i) for i in range(11)),
        )
        uns = {
            "title": "my test dataset",
            "X_approximate_distribution": "normal",
            "batch_condition": np.array({"batchA", "batchB"}),
            "schema_version": "3.0.0",
        }

        var = pandas.DataFrame(
            data=["gene", "spike-in", "gene"], columns=["feature_biotype"], index=non_zeros_X_layer_df.columns
        )

        adata = anndata.AnnData(
            X=non_zeros_X_layer_df, obs=obs, uns=uns, var=var, layers={"my_awesome_wonky_layer": zeros_layer_df}
        )
        adata_raw = anndata.AnnData(X=zeros_layer_df, obs=obs, uns=uns)
        adata.raw = adata_raw

        mock_read_h5ad.return_value = adata

        # Run the extraction method
        extracted_metadata = extract_metadata("dummy")

        # Verify that the "my_awesome_wonky_layer" was read and not the default X layer. The layer contains only zeros
        # which should result in a mean_genes_per_cell value of 0 compared to 3 if the X layer was read.
        self.assertEqual(extracted_metadata["mean_genes_per_cell"], 0)

    def mock_downloader_function(self, url, local_path, tracker, chunk_size):
        time.sleep(1)
        dataset = Dataset.get(self.session, self.dataset_id)
        dataset.update(tombstone=True)
        for x in range(10):
            if tracker.stop_downloader.is_set():
                return
            time.sleep(3)

    @patch("backend.dataset_pipeline.processing.download.downloader")
    @patch("backend.dataset_pipeline.processing.process_download_validate.from_url")
    def test__dataset_tombstoned_while_uploading(self, mock_from_url, mock_downloader):
        class file_url:
            scheme = "https"
            url = "url.com"

            @classmethod
            def file_info(cls):
                return {"size": 12}

        mock_downloader.side_effect = self.mock_downloader_function
        mock_from_url.return_value = file_url
        self.dataset_id = self.generate_dataset(self.session).id
        start = time.time()
        # check that changing the db status leads to an exception being raised
        with self.assertRaises(ProcessingCancelled):
            download_from_source_uri(
                self.dataset_id,
                "dropbox.com",
                "raw.h5ad",
            )
        end = time.time()
        # check that tombstoning ends the download thread early
        self.assertLess(end - start, 11)

    def test__download_from_source_uri_with_unhandled_scheme__raises_error(self):
        unhandled_uri = "unhandled_scheme://blah/foo"
        self.dataset_id = self.generate_dataset(self.session).id

        with self.assertRaises(ValueError):
            download_from_source_uri(self.dataset_id, unhandled_uri, "raw.h5ad")

    @mock_s3
    @patch("backend.dataset_pipeline.processing.process_download_validate.download_from_s3")
    def test__download_from_source_uri_with_s3_scheme__downloads_from_s3(self, mock_download_from_s3):
        test_dataset_id = self.generate_dataset(self.session).id

        bucket = "bucket"
        key = "key"
        local_file = "local.h5ad"
        s3_uri = f"s3://{bucket}/{key}"

        download_from_source_uri(test_dataset_id, s3_uri, local_file)

        mock_download_from_s3.assert_called_with(bucket_name=bucket, object_key=key, local_filename=local_file)
        dataset = Dataset.get(self.session, test_dataset_id)
        self.assertEqual(UploadStatus.UPLOADED, dataset.processing_status.upload_status)
