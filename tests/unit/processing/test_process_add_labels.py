import tempfile
from unittest.mock import MagicMock, Mock, patch

import anndata

from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_add_labels import ProcessAddLabels
from tests.unit.processing.base_processing_test import BaseProcessingTest


class ProcessingTest(BaseProcessingTest):
    @patch("scanpy.read_h5ad")
    def test_process_add_labels(self, mock_read_h5ad):
        """
        ProcessDownloadValidate should:
        1. Download the h5ad artifact
        2. Add labels to h5ad
        2. Set upload status to UPLOADED
        3. set h5ad status to UPLOADED
        4. upload the labeled file to S3
        """
        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0"
        self.crossref_provider.fetch_metadata = Mock(return_value=({}, "12.2345", 17169328.664))

        collection = self.generate_unpublished_collection(
            links=[Link(name=None, type="DOI", uri="http://doi.org/12.2345")]
        )
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        # This is where we're at when we start the SFN

        mock_read_h5ad.return_value = MagicMock(uns=dict())

        # TODO: ideally use a real h5ad
        processor = ProcessAddLabels(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        processor.extract_metadata = Mock()
        processor.populate_dataset_citation = Mock()
        processor.process(collection.version_id, dataset_version_id, "fake_bucket_name", "fake_datasets_bucket")

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)

        # Verify that the labeled (local.h5ad) file is there
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))
        # Verify that the labeled file is uploaded to the datasets bucket
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.h5ad"))

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(1, len(artifacts))
        artifact = artifacts[0]
        artifact.type = DatasetArtifactType.H5AD
        artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"

    def test_populate_dataset_citation__with_publication_doi(self):
        mock_adata = anndata.AnnData(X=None, obs=None, obsm=None, uns={}, var=None)
        self.crossref_provider.fetch_metadata = Mock(return_value=({}, "12.2345", 17169328.664))
        collection = self.generate_unpublished_collection(
            links=[Link(name=None, type="DOI", uri="https://doi.org/12.2345")]
        )
        with tempfile.NamedTemporaryFile(suffix=".h5ad") as f:
            mock_adata.write_h5ad(f.name)
            pal = ProcessAddLabels(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            dataset_version_id = DatasetVersionId()
            pal.populate_dataset_citation(collection.version_id, dataset_version_id, f.name)
            citation_str = (
                f"Publication: https://doi.org/12.2345 "
                f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
                f"CZ CELLxGENE Discover in Collection: https://domain/collections/{collection.collection_id}"
            )
            adata = anndata.read_h5ad(f.name)
            self.assertEqual(adata.uns["citation"], citation_str)

    def test_populate_dataset_citation__no_publication_doi(self):
        mock_adata = anndata.AnnData(X=None, obs=None, obsm=None, uns={}, var=None)
        collection = self.generate_unpublished_collection()
        with tempfile.NamedTemporaryFile(suffix=".h5ad") as f:
            mock_adata.write_h5ad(f.name)
            pal = ProcessAddLabels(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            dataset_version_id = DatasetVersionId()
            pal.populate_dataset_citation(collection.version_id, dataset_version_id, f.name)
            citation_str = (
                f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
                f"CZ CELLxGENE Discover in Collection: https://domain/collections/{collection.collection_id}"
            )
            adata = anndata.read_h5ad(f.name)
            self.assertEqual(adata.uns["citation"], citation_str)

    def test_process_add_labels_fail(self):
        """
        If the validation is not successful, the processing pipeline should:
        1. Set the processing status to INVALID
        2. Set a validation message accordingly
        """
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        manifest = IngestionManifest(anndata=dropbox_uri)

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        self.schema_validator.add_labels = Mock(side_effect=ValueError("Add labels error"))
        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pm.process(
            collection.version_id,
            dataset_version_id,
            "add_labels",
            manifest,
            "fake_bucket_name",
            "fake_datasets_bucket",
            "fake_cxg_bucket",
        )

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.FAILED)
