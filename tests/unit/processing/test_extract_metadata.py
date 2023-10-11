from unittest.mock import patch

import anndata
import numpy as np
import pandas

from backend.layers.common.entities import OntologyTermId, TissueOntologyTermId
from backend.layers.processing.process_download_validate import ProcessDownloadValidate
from tests.unit.processing.base_processing_test import BaseProcessingTest


class TestProcessingDownloadValidate(BaseProcessingTest):
    def setUp(self):
        super().setUp()
        self.pdv = ProcessDownloadValidate(
            self.business_logic, self.uri_provider, self.s3_provider, self.downloader, self.schema_validator
        )

        def mock_config_fn(name):
            if name == "schema_4_feature_flag":
                return "True"

        mock_config = patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
        mock_config.start()

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
                    np.array([["organoid", "tissue"][i] for i in tissue]).reshape(50001, 1),
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
                "tissue_type",
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

        extracted_metadata = self.pdv.extract_metadata("dummy")

        self.assertEqual(extracted_metadata.organism, [OntologyTermId("Homo sapiens", "NCBITaxon:8505")])

        self.assertCountEqual(
            extracted_metadata.tissue,
            [
                TissueOntologyTermId("lung", "UBERON:01", "organoid"),
                TissueOntologyTermId("liver", "UBERON:10", "tissue"),
            ],
        )

        self.assertCountEqual(
            extracted_metadata.assay,
            [
                OntologyTermId("10x", "EFO:001"),
                OntologyTermId("smartseq", "EFO:010"),
                OntologyTermId("cite-seq", "EFO:011"),
            ],
        )

        self.assertCountEqual(extracted_metadata.disease, [OntologyTermId("healthy", "MONDO:123")])

        self.assertCountEqual(
            extracted_metadata.sex,
            [OntologyTermId("male", "M"), OntologyTermId("female", "F"), OntologyTermId("fixed", "MF")],
        )

        self.assertCountEqual(
            extracted_metadata.self_reported_ethnicity,
            [OntologyTermId("solomon islander", "HANCESTRO:321"), OntologyTermId("orcadian", "HANCESTRO:456")],
        )

        self.assertCountEqual(
            extracted_metadata.development_stage,
            [
                OntologyTermId("adult", "HsapDv:0"),
                OntologyTermId("baby", "HsapDv:1"),
                OntologyTermId("tween", "HsapDv:2"),
            ],
        )

        self.assertCountEqual(
            extracted_metadata.suspension_type,
            ["cell", "nucleus", "na"],
        )

        self.assertCountEqual(
            extracted_metadata.donor_id,
            ["F1", "F2"],
        )

        self.assertEqual(extracted_metadata.x_approximate_distribution, "NORMAL")
        self.assertEqual(extracted_metadata.batch_condition, np.array({"batchA", "batchB"}))
        self.assertEqual(extracted_metadata.schema_version, "3.0.0")

        self.assertEqual(extracted_metadata.cell_count, 50001)
        self.assertEqual(extracted_metadata.primary_cell_count, 0)

        filter = np.where(adata.var.feature_biotype == "gene")[0]
        self.assertAlmostEqual(extracted_metadata.mean_genes_per_cell, np.count_nonzero(adata.X[:, filter]) / 50001)

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
                    np.array([["organoid", "tissue"][i] for i in tissue]).reshape(11, 1),
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
                "tissue_type",
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
        extracted_metadata = self.pdv.extract_metadata("dummy")

        # Verify that the "my_awesome_wonky_layer" was read and not the default X layer. The layer contains only zeros
        # which should result in a mean_genes_per_cell value of 0 compared to 3 if the X layer was read.
        self.assertEqual(extracted_metadata.mean_genes_per_cell, 0)
