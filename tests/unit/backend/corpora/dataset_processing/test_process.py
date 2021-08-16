import enum

import logging
import os
import pathlib
import shutil
import tempfile
import time
from unittest.mock import patch

import anndata
import boto3
import numpy
import pandas
from moto import mock_s3

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    DatasetArtifactType,
    DatasetArtifactFileType,
    UploadStatus,
    ValidationStatus,
    ConversionStatus,
)
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled
from backend.corpora.dataset_processing import process
from backend.corpora.dataset_processing.process import (
    convert_file_ignore_exceptions,
    download_from_dropbox_url,
)
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDatasetProcessing(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_dir = tempfile.mkdtemp()
        cls.h5ad_filename = pathlib.Path(cls.tmp_dir, "test.h5ad")
        cls.seurat_filename = pathlib.Path(cls.tmp_dir, "test.rds")
        cls.loom_filename = pathlib.Path(cls.tmp_dir, "test.loom")
        cls.cxg_filename = pathlib.Path(cls.tmp_dir, "test.cxg")

        cls.h5ad_filename.touch()
        cls.seurat_filename.touch()
        cls.loom_filename.touch()
        cls.cxg_filename.touch()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        shutil.rmtree(cls.tmp_dir)

    def setup_s3_bucket(self, bucket_name):
        # Mock S3 service if we don't have a mock api already running
        if os.getenv("BOTO_ENDPOINT_URL"):
            s3_args = {"endpoint_url": os.getenv("BOTO_ENDPOINT_URL")}
        else:
            s3_mock = mock_s3()
            s3_mock.start()
            s3_args = {}
            self.addCleanup(s3_mock.stop)
        s3 = boto3.client("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)
        self.s3_resource = boto3.resource("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)
        try:
            s3.create_bucket(
                Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
            )
        except self.s3_resource.meta.client.exceptions.BucketAlreadyExists:
            pass
        return s3

    def delete_s3_bucket(self, bucket_name):
        bucket = self.s3_resource.Bucket(bucket_name)
        if bucket.creation_date is not None:
            bucket.objects.all().delete()
            bucket.delete()

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
        sx = dev = numpy.random.choice([0, 1, 2], size=(50001))

        obs = pandas.DataFrame(
            numpy.hstack(
                [
                    numpy.array([["lung", "liver"][i] for i in tissue]).reshape(50001, 1),
                    numpy.array([["UBERON:01", "UBERON:10"][i] for i in tissue]).reshape(50001, 1),
                    numpy.array([["10x", "smartseq", "cite-seq"][i] for i in assay]).reshape(50001, 1),
                    numpy.array([["EFO:001", "EFO:010", "EFO:011"][i] for i in assay]).reshape(50001, 1),
                    numpy.random.choice(["healthy"], size=(50001, 1)),
                    numpy.random.choice(["MONDO:123"], size=(50001, 1)),
                    numpy.array([["male", "female", "fixed"][i] for i in sx]).reshape(50001, 1),
                    numpy.array([["M", "F", "MF"][i] for i in sx]).reshape(50001, 1),
                    numpy.array([["solomon islander", "orcadian"][i] for i in eth]).reshape(50001, 1),
                    numpy.array([["HANCESTRO:321", "HANCESTRO:456"][i] for i in eth]).reshape(50001, 1),
                    numpy.array([["adult", "baby", "tween"][i] for i in dev]).reshape(50001, 1),
                    numpy.array([["HsapDv:0", "HsapDv:1", "HsapDv:2"][i] for i in dev]).reshape(50001, 1),
                    numpy.random.choice(["Homo sapiens"], size=(50001, 1)),
                    numpy.random.choice(["NCBITaxon:8505"], size=(50001, 1)),
                    numpy.random.choice([0], size=(50001, 1)),
                    numpy.random.choice(["liver"], size=(50001, 1)),
                    numpy.random.choice(["Hepatic-1A"], size=(50001, 1)),
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
                "ethnicity",
                "ethnicity_ontology_term_id",
                "development_stage",
                "development_stage_ontology_term_id",
                "organism",
                "organism_ontology_term_id",
                "is_primary_data",
                "cell_type",
                "cell_type_ontology_term_id"
            ],
            index=(str(i) for i in range(50001)),
        )
        uns = {
            "title": "my test dataset",
            "X_normalization": "normal",
            "X_approximate_distribution": "normal",
        }

        adata = anndata.AnnData(X=df, obs=obs, uns=uns)
        mock_read_h5ad.return_value = adata

        extracted_metadata = process.extract_metadata("dummy")
        lab, ont = "label", "ontology_term_id"

        def list_equal(list1, list2, cmp_func):
            self.assertEqual(len(list1), len(list2))
            for el1 in list1:
                self.assertIn(el1, list2)
                el2 = list2[list2.index(el1)]
                if cmp_func:
                    cmp_func(el1, el2)

        list_equal(
            extracted_metadata["organism"], 
            [{lab: "Homo sapiens", ont: "NCBITaxon:8505"}], 
            self.assertDictEqual
        )

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

    @patch("scanpy.read_h5ad")
    def test_extract_metadata_find_raw_layer(self, mock_read_h5ad):
        # Setup anndata to be read
        non_zeros_X_layer_df = pandas.DataFrame(
            numpy.full((11, 3), 2), columns=list("ABC"), index=(str(i) for i in range(11))
        )
        zeros_layer_df = pandas.DataFrame(numpy.zeros((11, 3)), columns=list("ABC"), index=(str(i) for i in range(11)))

        obs = pandas.DataFrame(
            numpy.hstack(
                [
                    numpy.array([["lung", "liver"][i] for i in numpy.random.choice([0, 1], size=(11))]).reshape(11, 1),
                    numpy.array(
                        [["UBERON:01", "UBERON:10"][i] for i in numpy.random.choice([0, 1], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["10x", "smartseq", "cite-seq"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["EFO:001", "EFO:010", "EFO:011"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.random.choice(["healthy"], size=(11, 1)),
                    numpy.random.choice(["MONDO:123"], size=(11, 1)),
                    numpy.array(
                        [["male", "female", "fixed"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["M", "F", "MF"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["solomon islander", "orcadian"][i] for i in numpy.random.choice([0, 1], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["HANCESTRO:321", "HANCESTRO:456"][i] for i in numpy.random.choice([0, 1], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["adult", "baby", "tween"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.array(
                        [["HsapDv:0", "HsapDv:1", "HsapDv:2"][i] for i in numpy.random.choice([0, 1, 2], size=(11))]
                    ).reshape(11, 1),
                    numpy.random.choice(["Homo sapiens"], size=(11, 1)),
                    numpy.random.choice(["NCBITaxon:8505"], size=(11, 1)),
                    numpy.random.choice([0], size=(11, 1)),
                    numpy.random.choice(["liver"], size=(11, 1)),
                    numpy.random.choice(["Hepatic-1A"], size=(11, 1)),
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
                "ethnicity",
                "ethnicity_ontology_term_id",
                "development_stage",
                "development_stage_ontology_term_id",
                "organism",
                "organism_ontology_term_id",
                "is_primary_data",
                "cell_type",
                "cell_type_ontology_term_id"
            ],
            index=(str(i) for i in range(11)),
        )
        uns = {
            "title": "my test dataset",
            "X_normalization": "normal",
            "X_approximate_distribution": "normal",
        }

        adata = anndata.AnnData(
            X=non_zeros_X_layer_df, obs=obs, uns=uns, layers={"my_awesome_wonky_layer": zeros_layer_df}
        )
        adata_raw = anndata.AnnData(
            X=zeros_layer_df, obs=obs, uns=uns
        )
        adata.raw = adata_raw
        
        mock_read_h5ad.return_value = adata

        # Run the extraction method
        extracted_metadata = process.extract_metadata("dummy")

        # Verify that the "my_awesome_wonky_layer" was read and not the default X layer. The layer contains only zeros
        # which should result in a mean_genes_per_cell value of 0 compared to 3 if the X layer was read.
        self.assertEqual(extracted_metadata["mean_genes_per_cell"], 0)

    def test_update_db(self):

        collection = Collection.create(self.session, visibility=CollectionVisibility.PRIVATE)
        dataset = Dataset.create(
            self.session, collection_id=collection.id, collection_visibility=CollectionVisibility.PRIVATE
        )
        dataset_id = dataset.id

        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        fake_env.start()

        process.update_db(
            dataset_id, metadata={"sex": ["male", "female"], "explorer_url": "https://cellxgene.com/data"}
        )
        self.session.expire(dataset)
        self.assertListEqual(Dataset.get(self.session, dataset_id).sex, ["male", "female"])

        artifact = {
            "filename": "test_filename",
            "filetype": DatasetArtifactFileType.H5AD,
            "type": DatasetArtifactType.REMIX,
            "user_submitted": True,
            "s3_uri": "s3://test_uri",
        }
        process.update_db(dataset_id, metadata={"artifacts": [artifact]})
        self.session.expire(dataset)
        self.assertEqual(len(Dataset.get(self.session, dataset_id).artifacts), 1)
        self.assertEqual(Dataset.get(self.session, dataset_id).artifacts[0].filename, "test_filename")
        self.assertEqual(Dataset.get(self.session, dataset_id).explorer_url, "https://cellxgene.com/data")

        process.update_db(
            dataset_id, processing_status={"upload_status": UploadStatus.UPLOADING, "upload_progress": 0.5}
        )
        self.session.expire(dataset)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_progress, 0.5)
        self.assertIsNone(Dataset.get(self.session, dataset_id).processing_status.validation_status)

        process.update_db(
            dataset_id,
            processing_status={
                "upload_status": UploadStatus.UPLOADED,
                "upload_progress": 1,
                "validation_status": ValidationStatus.VALIDATING,
            },
        )
        self.session.expire(dataset)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_status, UploadStatus.UPLOADED)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_progress, 1)
        self.assertEqual(
            Dataset.get(self.session, dataset_id).processing_status.validation_status, ValidationStatus.VALIDATING
        )

        fake_env.stop()

    def test_update_db__tombstoned_dataset(self):
        dataset = self.generate_dataset(self.session, tombstone=True)
        dataset_id = dataset.id

        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        fake_env.start()

        with self.assertRaises(ProcessingCancelled):
            process.update_db(dataset_id, metadata={"sex": ["male", "female"]})

    @patch("backend.corpora.dataset_processing.process.make_cxg")
    @patch("backend.corpora.dataset_processing.process.subprocess.run")
    def test_create_explorer_cxg(self, mock_subprocess, mock_cxg):
        mock_cxg.return_value = str(self.cxg_filename)
        dataset = self.generate_dataset(self.session)
        dataset_id = dataset.id
        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        fake_env.start()

        explorer_bucket = "CELLXGENE-HOSTED-TEST"

        process.process_cxg(str(self.h5ad_filename), dataset_id, explorer_bucket)

        dataset = Dataset.get(self.session, dataset_id)
        artifacts = dataset.artifacts

        self.assertEqual(len(artifacts), 1)
        self.assertEqual(artifacts[0].dataset_id, dataset_id)
        self.assertEqual(artifacts[0].s3_uri, f"s3://{explorer_bucket}/{dataset_id}.cxg/")
        self.assertEqual(artifacts[0].filetype, DatasetArtifactFileType.CXG)

    @patch("backend.corpora.dataset_processing.process.make_loom")
    @patch("backend.corpora.dataset_processing.process.make_seurat")
    def test_create_artifacts(self, make_seurat, make_loom):
        make_loom.return_value = str(self.loom_filename)
        make_seurat.return_value = str(self.seurat_filename)
        artifact_bucket = "test-artifact-bucket"
        test_dataset_id = self.generate_dataset(
            self.session,
        ).id

        s3 = self.setup_s3_bucket(artifact_bucket)

        bucket_prefix = process.get_bucket_prefix(test_dataset_id)
        process.create_artifacts(str(self.h5ad_filename), test_dataset_id, artifact_bucket)
        dataset = Dataset.get(self.session, test_dataset_id)
        artifacts = dataset.artifacts
        processing_status = dataset.processing_status

        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_loom_status)
        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_rds_status)
        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_anndata_status)

        self.assertEqual(len(artifacts), 3)

        self.assertTrue(all(a.user_submitted for a in artifacts))
        self.assertTrue(all(a.s3_uri.startswith(f"s3://{artifact_bucket}/{bucket_prefix}/") for a in artifacts))
        self.assertEqual(len(set(a.filetype for a in artifacts)), 3)
        self.assertTrue(all(a.type == DatasetArtifactType.REMIX for a in artifacts))

        resp = s3.list_objects_v2(Bucket=artifact_bucket, Prefix=bucket_prefix)
        s3_filenames = [os.path.basename(c["Key"]) for c in resp["Contents"]]
        self.assertEqual(len(s3_filenames), 3)
        self.assertIn(str(self.h5ad_filename.parts[-1]), s3_filenames)
        self.assertIn(str(self.seurat_filename.parts[-1]), s3_filenames)
        self.assertIn(str(self.loom_filename.parts[-1]), s3_filenames)

        # cleanup
        self.delete_s3_bucket(artifact_bucket)

    def test__create_artifact__negative(self):
        artifact_bucket = "test-artifact-bucket"
        test_dataset = self.generate_dataset(
            self.session,
        )
        bucket_prefix = process.get_bucket_prefix(test_dataset.id)
        self.setup_s3_bucket(artifact_bucket)

        with self.subTest("file does not exist"):
            self.assertRaises(
                TypeError,
                process.create_artifact,
                None,
                DatasetArtifactFileType.H5AD,
                bucket_prefix,
                test_dataset.id,
                artifact_bucket,
            )

        with self.subTest("invalid artifact type"):

            class BadEnum(enum.Enum):
                fake = "fake"

            self.assertRaises(
                CorporaException,
                process.create_artifact,
                str(self.h5ad_filename),
                BadEnum.fake,
                bucket_prefix,
                test_dataset.id,
                artifact_bucket,
            )

        with self.subTest("dataset does not exist"):
            self.assertRaises(
                CorporaException,
                process.create_artifact,
                str(self.h5ad_filename),
                DatasetArtifactFileType.H5AD,
                process.get_bucket_prefix("1234"),
                "1234",
                artifact_bucket,
            )

        with self.subTest("bucket does not exist"):
            self.assertRaises(
                boto3.exceptions.S3UploadFailedError,
                process.create_artifact,
                str(self.h5ad_filename),
                DatasetArtifactFileType.H5AD,
                bucket_prefix,
                test_dataset.id,
                "fake-bucket",
            )

        # cleanup
        self.delete_s3_bucket(artifact_bucket)

    @patch("backend.corpora.dataset_processing.process.make_loom")
    @patch("backend.corpora.dataset_processing.process.make_seurat")
    def test_process_continues_with_loom_conversion_failures(self, mock_seurat, mock_loom):
        mock_loom.side_effect = RuntimeError("Loom conversion failed")
        mock_seurat.return_value = str(self.h5ad_filename).replace(".h5ad", ".rds")
        test_dataset_id = self.generate_dataset(
            self.session,
        ).id
        bucket_prefix = process.get_bucket_prefix(test_dataset_id)
        artifact_bucket = "test-artifact-bucket"
        s3 = self.setup_s3_bucket(artifact_bucket)
        process.create_artifacts(str(self.h5ad_filename), test_dataset_id, artifact_bucket)
        dataset = Dataset.get(self.session, test_dataset_id)
        artifacts = dataset.artifacts
        processing_status = dataset.processing_status

        self.assertEqual(ConversionStatus.FAILED, processing_status.conversion_loom_status)
        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_rds_status)
        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_anndata_status)

        self.assertEqual(len(artifacts), 2)
        resp = s3.list_objects_v2(Bucket=artifact_bucket, Prefix=bucket_prefix)
        s3_filenames = [os.path.basename(c["Key"]) for c in resp["Contents"]]
        self.assertEqual(len(s3_filenames), 2)
        self.assertNotIn(str(self.loom_filename.parts[-1]), s3_filenames)

        # cleanup
        self.delete_s3_bucket(artifact_bucket)

    @patch("backend.corpora.dataset_processing.process.make_loom")
    @patch("backend.corpora.dataset_processing.process.make_seurat")
    def test_process_continues_with_seurat_conversion_failures(self, mock_seurat, mock_loom):
        mock_loom.return_value = str(self.loom_filename)
        mock_seurat.side_effect = RuntimeError("seurat conversion failed")
        test_dataset_id = self.generate_dataset(
            self.session,
        ).id
        bucket_prefix = process.get_bucket_prefix(test_dataset_id)
        artifact_bucket = "test-artifact-bucket"
        s3 = self.setup_s3_bucket(artifact_bucket)
        process.create_artifacts(str(self.h5ad_filename), test_dataset_id, artifact_bucket)
        dataset = Dataset.get(self.session, test_dataset_id)
        artifacts = dataset.artifacts
        processing_status = dataset.processing_status

        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_loom_status)
        self.assertEqual(ConversionStatus.FAILED, processing_status.conversion_rds_status)
        self.assertEqual(ConversionStatus.CONVERTED, processing_status.conversion_anndata_status)

        self.assertEqual(len(artifacts), 2)
        resp = s3.list_objects_v2(Bucket=artifact_bucket, Prefix=bucket_prefix)
        s3_filenames = [os.path.basename(c["Key"]) for c in resp["Contents"]]
        self.assertEqual(len(s3_filenames), 2)
        self.assertNotIn(str(self.seurat_filename.parts[-1]), s3_filenames)

        # cleanup
        self.delete_s3_bucket(artifact_bucket)

    @patch("backend.corpora.dataset_processing.process.make_cxg")
    def test_process_continues_with_cxg_conversion_failures(self, mock_cxg):
        mock_cxg.side_effect = RuntimeError("cxg conversion failed")
        test_dataset_id = self.generate_dataset(
            self.session,
        ).id
        artifact_bucket = "test-artifact-bucket"
        process.process_cxg(str(self.h5ad_filename), test_dataset_id, artifact_bucket)
        dataset = Dataset.get(self.session, test_dataset_id)
        processing_status = dataset.processing_status
        self.assertEqual(ConversionStatus.FAILED, processing_status.conversion_cxg_status)

    def test__convert_file_ignore_exceptions__fail(self):
        def converter(_file):
            raise RuntimeError("conversion_failed")

        with self.assertLogs(process.logger, logging.ERROR):
            filename, status = convert_file_ignore_exceptions(converter, self.h5ad_filename, "error")
        self.assertIsNone(filename)
        self.assertEqual(ConversionStatus.FAILED, status)

    def mock_downloader_function(self, url, local_path, tracker, chunk_size):
        time.sleep(1)
        dataset = Dataset.get(self.session, self.dataset_id)
        dataset.update(tombstone=True)
        for x in range(10):
            if tracker.stop_downloader.is_set():
                return
            time.sleep(3)

    @patch("backend.corpora.dataset_processing.download.downloader")
    @patch("backend.corpora.dataset_processing.process.from_url")
    def test__dataset_tombstoned_while_uploading(self, mock_from_url, mock_downloader):
        class file_info:
            @classmethod
            def file_info(self):
                return {"size": 12}

            @property
            def url(self):
                return "url.com"

        mock_downloader.side_effect = self.mock_downloader_function
        mock_from_url.return_value = file_info
        self.dataset_id = self.generate_dataset(self.session).id
        start = time.time()
        # check that changing the db status leads to an exception being raised
        with self.assertRaises(ProcessingCancelled):
            download_from_dropbox_url(
                self.dataset_id,
                "dropbox.com",
                "local.h5ad",
            )
        end = time.time()
        # check that tombstoning ends the download thread early
        self.assertLess(end - start, 11)

    @patch("backend.corpora.dataset_processing.process.make_cxg")
    @patch("backend.corpora.dataset_processing.process.download_from_dropbox_url")
    @patch("backend.corpora.dataset_processing.process.validate_h5ad_file")
    @patch("backend.corpora.dataset_processing.process.extract_metadata")
    def test__cxg_not_created_when_metadata_extraction_fails(
        self, mock_metadata_extraction, mock_validate_h5ad, mock_dropbox_download, mock_cxg
    ):
        mock_metadata_extraction.side_effect = RuntimeError("metadata extraction failed")
        mock_dropbox_download.return_value = self.h5ad_filename
        mock_cxg.return_value = str(self.cxg_filename)
        dataset = self.generate_dataset(self.session)
        dataset_id = dataset.id

        explorer_bucket = "CELLXGENE-HOSTED-TEST"
        artifact_bucket = "test-artifact-bucket"

        with self.assertRaises(RuntimeError):
            process.process(dataset_id, str(self.h5ad_filename), explorer_bucket, artifact_bucket)

        dataset = Dataset.get(self.session, dataset_id)
        processing_status = dataset.processing_status
        self.assertEqual(None, processing_status.conversion_cxg_status)
