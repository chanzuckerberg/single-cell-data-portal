import os
import unittest
from unittest.mock import patch

import anndata
import numpy
import pandas

from backend.corpora.common.corpora_orm import CollectionVisibility, DatasetArtifactType, DatasetArtifactFileType
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.entities.dataset import Dataset

from backend.corpora.dataset_processing import process


class TestDatasetProcessing(unittest.TestCase):
    @patch.dict(
        os.environ,
        {
            "DROPBOX_URL": "xxx",
            "ARTIFACT_BUCKET": "yyy",
            "CELLXGENE_BUCKET": "zzz",
            "DATASET_ID": "aaa",
            "DEPLOYMENT_STAGE": "test",
        },
    )
    def test_check_env_good(self):
        process.check_env()

    @patch.dict(os.environ, {"DROPBOX_URL": "xxx", "ARTIFACT_BUCKET": "yyy", "CELLXGENE_BUCKET": "zzz"})
    def test_check_env_bad(self):
        with self.assertRaises(EnvironmentError):
            process.check_env()

    @patch("scanpy.read_h5ad")
    def test_extract_metadata(self, mock_read_h5ad):

        df = pandas.DataFrame(
            numpy.random.randint(10, size=(50001, 5)) * 50, columns=list("ABCDE"), index=(str(i) for i in range(50001))
        )

        tissue = numpy.random.choice([0, 1], size=(50001))
        assay = numpy.random.choice([0, 1, 2], size=(50001))
        eth = numpy.random.choice([0, 1], size=(50001))
        dev = numpy.random.choice([0, 1, 2], size=(50001))

        obs = pandas.DataFrame(
            numpy.hstack(
                [
                    numpy.array([["lung", "liver"][i] for i in tissue]).reshape(50001, 1),
                    numpy.array([["UBERON:01", "UBERON:10"][i] for i in tissue]).reshape(50001, 1),
                    numpy.array([["10x", "smartseq", "cite-seq"][i] for i in assay]).reshape(50001, 1),
                    numpy.array([["EFO:001", "EFO:010", "EFO:011"][i] for i in assay]).reshape(50001, 1),
                    numpy.random.choice(["healthy"], size=(50001, 1)),
                    numpy.random.choice(["MONDO:123"], size=(50001, 1)),
                    numpy.random.choice(["male", "female"], size=(50001, 1)),
                    numpy.array([["solomon islander", "orcadian"][i] for i in eth]).reshape(50001, 1),
                    numpy.array([["HANCESTRO:321", "HANCESTRO:456"][i] for i in eth]).reshape(50001, 1),
                    numpy.array([["adult", "baby", "tween"][i] for i in dev]).reshape(50001, 1),
                    numpy.array([["HsapDv:0", "HsapDv:1", "HsapDv:2"][i] for i in dev]).reshape(50001, 1),
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
                "ethnicity",
                "ethnicity_ontology_term_id",
                "development_stage",
                "development_stage_ontology_term_id",
            ],
            index=(str(i) for i in range(50001)),
        )
        uns = {
            "organism": "Homo sapiens",
            "organism_ontology_term_id": "NCBITaxon:8505",
            "layer_descriptions": {"X": "raw"},
        }

        adata = anndata.AnnData(X=df, obs=obs, uns=uns)
        mock_read_h5ad.return_value = adata

        extracted_metadata = process.extract_metadata("dummy")
        lab, ont = "label", "ontology_term_id"

        self.assertDictEqual(extracted_metadata["organism"], {lab: "Homo sapiens", ont: "NCBITaxon:8505"})

        def list_equal(list1, list2, cmp_func):
            self.assertEqual(len(list1), len(list2))
            for el1 in list1:
                self.assertIn(el1, list2)
                el2 = list2[list2.index(el1)]
                if cmp_func:
                    cmp_func(el1, el2)

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

        self.assertListEqual(sorted(extracted_metadata["sex"]), sorted(["male", "female"]))

        list_equal(
            extracted_metadata["ethnicity"],
            [{lab: "solomon islander", ont: "HANCESTRO:321"}, {lab: "orcadian", ont: "HANCESTRO:456"}],
            self.assertDictEqual,
        )

        list_equal(
            extracted_metadata["development_stage"],
            [{lab: "adult", ont: "HsapDv:0"}, {lab: "baby", ont: "HsapDv:1"}, {lab: "tween", ont: "HsapDv:2"}],
            self.assertDictEqual,
        )
        self.assertEqual(extracted_metadata["cell_count"], 50001)
        self.assertAlmostEqual(extracted_metadata["mean_genes_per_cell"], numpy.count_nonzero(df) / 50001)

    def test_update_db(self):

        collection = Collection.create(visibility=CollectionVisibility.PUBLIC)
        dataset = Dataset.create(collection_id=collection.id, collection_visibility=CollectionVisibility.PUBLIC)
        dataset_id = dataset.id

        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        Dataset.db.session.expire_all()
        fake_env.start()

        process.update_db(metadata={"sex": ["male", "female"]})
        Dataset.db.session.expire_all()

        self.assertListEqual(Dataset.get(dataset_id).sex, ["male", "female"])

        artifact = {
            "filename": "test_filename",
            "filetype": DatasetArtifactFileType.H5AD,
            "type": DatasetArtifactType.REMIX,
            "user_submitted": True,
            "s3_uri": "s3://test_uri",
        }
        dep_dir = {"url": "https://cellxgene.com/data"}
        process.update_db(metadata={"artifacts": [artifact], "deployment_directories": [dep_dir]})
        Dataset.db.session.expire_all()

        self.assertEqual(len(Dataset.get(dataset_id).artifacts), 1)
        self.assertEqual(Dataset.get(dataset_id).artifacts[0].filename, "test_filename")

        self.assertEqual(Dataset.get(dataset_id).deployment_directories[0].url, "https://cellxgene.com/data")

        fake_env.stop()


class TestDropBox(unittest.TestCase):
    def test_fix_url_direct_download(self):

        self.assertEqual(
            process.fix_dropbox_url("https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1"),
            "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
        )

        self.assertEqual(
            process.fix_dropbox_url("https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx"),
            "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
        )

        self.assertEqual(
            process.fix_dropbox_url("https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=0"),
            "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
        )

        self.assertEqual(
            process.fix_dropbox_url("https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever"),
            "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
        )

    def test_fix_url_bad(self):

        self.assertIsNone(process.fix_dropbox_url("https://www.googledrive.com/s/a1b/example.docx"))
        self.assertIsNone(
            process.fix_dropbox_url("http://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1")
        )

    @patch("subprocess.run")
    def test_fetch_dropbox_url_ok(self, mock_run):

        full_dp_url = "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1"
        process.fetch_dropbox_url(full_dp_url, "my_file.h5ad")

        print(mock_run.call_args[0])

        self.assertEqual(mock_run.call_args[0], (["wget", "-nv", full_dp_url, "-O", "my_file.h5ad"],))
