import json
import os
import tempfile
from unittest.mock import Mock, patch

import anndata
import pytest
import scanpy
import tiledb
from parameterized import parameterized

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.processing.dataset_metadata_update import DatasetMetadataUpdater
from backend.layers.processing.exceptions import ProcessingFailed
from backend.layers.processing.utils.cxg_generation_utils import convert_dictionary_to_cxg_group
from backend.layers.thirdparty.s3_provider_mock import MockS3Provider
from tests.unit.backend.layers.common.base_test import DatasetArtifactUpdate, DatasetStatusUpdate
from tests.unit.processing.base_processing_test import BaseProcessingTest


class TestUpdateMetadataHandler(BaseProcessingTest):
    def setUp(self):
        super().setUp()
        self.business_logic.s3_provider = MockS3Provider()
        self.updater = DatasetMetadataUpdater(
            self.business_logic, "artifact_bucket", "cellxgene_bucket", "datasets_bucket", "spatial_deep_zoom_dir"
        )
        self.updater.update_raw_h5ad = Mock()
        self.updater.update_h5ad = Mock()
        self.updater.update_cxg = Mock()

        def mock_download(source_uri, local_path):
            return local_path

        self.updater.download_from_source_uri = Mock(side_effect=mock_download)

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata(self, *args):
        current_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )
        self.updater.has_valid_artifact_statuses = Mock(return_value=True)
        self.updater.update_metadata(
            current_dataset_version_id, new_dataset_version_id, DatasetArtifactMetadataUpdate(citation="New Citation")
        )

        # skip raw_h5ad update since no updated fields are expected fields in raw H5AD
        self.updater.update_raw_h5ad.assert_not_called()
        self.updater.update_h5ad.assert_called_once()
        self.updater.update_cxg.assert_called_once()

        # check that collection version maps to dataset version with updated metadata
        collection_version = self.business_logic.get_collection_version(collection_version_id)
        new_dataset_version = collection_version.datasets[0]
        new_dataset_version_id = new_dataset_version.version_id
        artifacts = [(artifact.uri, artifact.type) for artifact in new_dataset_version.artifacts]
        assert (f"s3://artifact_bucket/{new_dataset_version_id}/raw.h5ad", DatasetArtifactType.RAW_H5AD) in artifacts

        assert new_dataset_version.status.upload_status == DatasetUploadStatus.UPLOADED
        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.SUCCESS

        # RDS should be skipped
        assert new_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED

        assert self.updater.s3_provider.uri_exists(f"s3://artifact_bucket/{new_dataset_version_id}/raw.h5ad")

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__rds_skipped(self, *args):
        current_dataset_version = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.RAW_H5AD, "s3://fake.bucket/raw.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "s3://fake.bucket/local.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "s3://fake.bucket/local.cxg"),
            ],
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.SKIPPED),
            ],
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )
        self.updater.has_valid_artifact_statuses = Mock(return_value=True)
        self.updater.update_metadata(
            current_dataset_version_id, new_dataset_version_id, DatasetArtifactMetadataUpdate(citation="New Citation")
        )

        self.updater.update_raw_h5ad.assert_not_called()
        self.updater.update_h5ad.assert_called_once()
        self.updater.update_cxg.assert_called_once()

        # check that collection version maps to dataset version with updated metadata
        collection_version = self.business_logic.get_collection_version(collection_version_id)
        new_dataset_version = collection_version.datasets[0]

        assert new_dataset_version.status.upload_status == DatasetUploadStatus.UPLOADED
        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.SUCCESS

        # RDS should be skipped
        assert new_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__raw_h5ad_updated(self, *args):
        current_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )
        self.updater.has_valid_artifact_statuses = Mock(return_value=True)
        self.updater.update_metadata(
            current_dataset_version_id, new_dataset_version_id, DatasetArtifactMetadataUpdate(title="New Dataset Title")
        )

        self.updater.update_raw_h5ad.assert_called_once()
        self.updater.update_h5ad.assert_called_once()
        self.updater.update_cxg.assert_called_once()

        # check that collection version maps to dataset version with updated metadata
        collection_version = self.business_logic.get_collection_version(collection_version_id)
        new_dataset_version = collection_version.datasets[0]

        # RDS should be skipped
        assert new_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED

        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.SUCCESS

    def test_update_metadata__current_dataset_version_bad_processing_status(self):
        current_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.FAILURE),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        self.updater.upload_raw_h5ad = Mock()
        self.updater.update_metadata(current_dataset_version_id, None, None)

        self.updater.upload_raw_h5ad.assert_not_called()

    def test_update_metadata__new_dataset_version_bad_processing_status(self):
        current_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        new_dataset_version_id = DatasetVersionId(new_dataset_version.dataset_version_id)
        self.updater.upload_raw_h5ad = Mock()
        self.updater.update_metadata(current_dataset_version_id, new_dataset_version_id, None)

        self.updater.upload_raw_h5ad.assert_not_called()

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__error_if_missing_raw_h5ad(self, *args):
        current_dataset_version = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "s3://fake.bucket/local.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "s3://fake.bucket/local.cxg"),
                DatasetArtifactUpdate(DatasetArtifactType.RDS, "s3://fake.bucket/local.rds"),
            ],
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ],
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )

        with pytest.raises(ValueError):
            self.updater.update_metadata(
                current_dataset_version_id,
                new_dataset_version_id,
                None,
            )

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__missing_labeled_h5ad(self, *args):
        current_dataset_version = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.RAW_H5AD, "s3://fake.bucket/raw.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "s3://fake.bucket/local.cxg"),
                DatasetArtifactUpdate(DatasetArtifactType.RDS, "s3://fake.bucket/local.rds"),
            ],
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ],
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )

        with pytest.raises(ProcessingFailed):
            self.updater.update_metadata(current_dataset_version_id, new_dataset_version_id, None)

        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)

        assert new_dataset_version.status.h5ad_status == DatasetConversionStatus.FAILED
        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.FAILURE

        # RDS should be skipped
        assert new_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__missing_cxg(self, *args):
        current_dataset_version = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.RAW_H5AD, "s3://fake.bucket/raw.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "s3://fake.bucket/local.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.RDS, "s3://fake.bucket/local.rds"),
            ],
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ],
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )

        with pytest.raises(ProcessingFailed):
            self.updater.update_metadata(current_dataset_version_id, new_dataset_version_id, None)

        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)

        assert new_dataset_version.status.cxg_status == DatasetConversionStatus.FAILED
        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.FAILURE

        # RDS should be skipped
        assert new_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED

    @patch("backend.common.utils.dl_sources.uri.downloader")
    def test_update_metadata__invalid_artifact_status(self, *args):
        current_dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
            ]
        )
        collection_version_id = CollectionVersionId(current_dataset_version.collection_version_id)
        current_dataset_version_id = DatasetVersionId(current_dataset_version.dataset_version_id)
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version_id,
            start_step_function=False,
        )
        self.updater.has_valid_artifact_statuses = Mock(return_value=False)
        with pytest.raises(ProcessingFailed):
            self.updater.update_metadata(current_dataset_version_id, new_dataset_version_id, None)

        collection_version = self.business_logic.get_collection_version(collection_version_id)
        new_dataset_version = collection_version.datasets[0]
        assert new_dataset_version.status.processing_status == DatasetProcessingStatus.FAILURE


class TestArtifactUpdates(BaseProcessingTest):
    def setUp(self):
        super().setUp()
        self.business_logic.s3_provider = MockS3Provider()
        self.updater = DatasetMetadataUpdater(
            self.business_logic, "artifact_bucket", "cellxgene_bucket", "datasets_bucket", "spatial_deep_zoom_dir"
        )
        mock_uns = {"title": "old title", "citation": "old citation", "schema_version": "3.0.0"}
        mock_adata = anndata.AnnData(X=None, obs=None, obsm=None, uns=mock_uns, var=None)
        self.temp_files = []

        def mock_download(source_uri, local_path):
            mock_adata.write_h5ad(local_path)
            self.temp_files.append(local_path)
            return local_path

        self.updater.download_from_source_uri = Mock(side_effect=mock_download)

    def takeDown(self):
        super().takeDown()
        for file in self.temp_files:
            os.remove(file)

    @patch("backend.layers.processing.dataset_metadata_update.os.remove")
    def test_update_raw_h5ad(self, *args):
        collection_version = self.generate_unpublished_collection(add_datasets=1)
        current_dataset_version = collection_version.datasets[0]
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version.version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version.version_id,
            start_step_function=False,
        )
        key_prefix = new_dataset_version_id.id
        metadata_update = DatasetArtifactMetadataUpdate(
            citation="Publication DOI www.doi.org/567.8", title="New Dataset Title"
        )

        self.updater.update_raw_h5ad(None, key_prefix, new_dataset_version_id, metadata_update)

        local_filename = CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME
        # check mock_anndata object
        mock_anndata = scanpy.read_h5ad(local_filename)
        assert mock_anndata.uns["citation"] == "Publication DOI www.doi.org/567.8"
        assert mock_anndata.uns["title"] == "New Dataset Title"
        assert mock_anndata.uns["schema_version"] == "3.0.0"  # not changed from set-up
        # check s3 uris exist
        assert self.updater.s3_provider.uri_exists(f"s3://artifact_bucket/{new_dataset_version_id.id}/{local_filename}")
        # check DB DatasetVersion
        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)
        artifacts = [(artifact.uri, artifact.type) for artifact in new_dataset_version.artifacts]
        assert (
            f"s3://artifact_bucket/{new_dataset_version_id.id}/{local_filename}",
            DatasetArtifactType.RAW_H5AD,
        ) in artifacts
        # check processing status
        assert new_dataset_version.status.upload_status == DatasetUploadStatus.UPLOADED

    @patch("backend.layers.processing.dataset_metadata_update.os.remove")
    def test_update_h5ad(self, *args):
        collection_version = self.generate_unpublished_collection(add_datasets=1)
        current_dataset_version = collection_version.datasets[0]
        new_dataset_version_id, _ = self.business_logic.ingest_dataset(
            collection_version_id=collection_version.version_id,
            url=None,
            file_size=0,
            current_dataset_version_id=current_dataset_version.version_id,
            start_step_function=False,
        )
        key_prefix = new_dataset_version_id.id
        metadata_update = DatasetArtifactMetadataUpdate(
            citation="Publication DOI www.doi.org/567.8", title="New Dataset Title", schema_reference="new field"
        )

        # Mock anndata object
        self.updater.update_h5ad(None, current_dataset_version, key_prefix, new_dataset_version_id, metadata_update)

        local_filename = CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME
        mock_anndata = scanpy.read_h5ad(local_filename)
        # check mock_anndata object
        assert mock_anndata.uns["citation"] == "Publication DOI www.doi.org/567.8"
        assert mock_anndata.uns["title"] == "New Dataset Title"
        assert mock_anndata.uns["schema_version"] == "3.0.0"
        assert mock_anndata.uns["schema_reference"] == "new field"
        # check s3 uris exist
        assert self.updater.s3_provider.uri_exists(f"s3://artifact_bucket/{new_dataset_version_id.id}/{local_filename}")
        assert self.updater.s3_provider.uri_exists(f"s3://datasets_bucket/{new_dataset_version_id.id}.h5ad")
        # check DB DatasetVersion
        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)
        assert new_dataset_version.metadata.citation == "Publication DOI www.doi.org/567.8"
        assert new_dataset_version.metadata.name == "New Dataset Title"
        assert new_dataset_version.metadata.schema_version == collection_version.datasets[0].metadata.schema_version
        artifacts = [(artifact.uri, artifact.type) for artifact in new_dataset_version.artifacts]
        assert (
            f"s3://artifact_bucket/{new_dataset_version_id.id}/{local_filename}",
            DatasetArtifactType.H5AD,
        ) in artifacts
        # check processing status
        assert new_dataset_version.status.validation_status == DatasetValidationStatus.VALID
        assert new_dataset_version.status.h5ad_status == DatasetConversionStatus.CONVERTED

    def test_update_cxg(self):
        with tempfile.TemporaryDirectory() as tempdir:
            testing_cxg_temp_directory = tempdir

            cxg_metadata = {
                "cxg_version": "1.0.0",
                "corpora": json.dumps(
                    {"schema_version": "3.0.0", "citation": "Publication www.doi.org/123.4", "title": "Dataset Title 1"}
                ),
            }
            convert_dictionary_to_cxg_group(
                testing_cxg_temp_directory, cxg_metadata, group_metadata_name="cxg_group_metadata"
            )

            collection_version = self.generate_unpublished_collection(add_datasets=1)
            current_dataset_version = collection_version.datasets[0]
            new_dataset_version_id, _ = self.business_logic.ingest_dataset(
                collection_version_id=collection_version.version_id,
                url=None,
                file_size=0,
                current_dataset_version_id=current_dataset_version.version_id,
                start_step_function=False,
            )
            metadata_update = DatasetArtifactMetadataUpdate(
                citation="Publication www.doi.org/567.8", title="New Dataset Title"
            )
            self.updater.update_cxg(
                None,
                testing_cxg_temp_directory,
                current_dataset_version.version_id,
                new_dataset_version_id,
                metadata_update,
            )

            # check new cxg directory exists
            assert self.updater.s3_provider.uri_exists(testing_cxg_temp_directory)
            # check spatial zoom directory was not created for non-spatial dataset
            assert not self.updater.s3_provider.uri_exists(
                f"s3://{self.updater.spatial_deep_zoom_dir}/{new_dataset_version_id.id}"
            )

            # check DB artifacts + status are updated
            new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)
            artifacts = [(artifact.uri, artifact.type) for artifact in new_dataset_version.artifacts]
            assert (testing_cxg_temp_directory, DatasetArtifactType.CXG) in artifacts

            assert new_dataset_version.status.cxg_status == DatasetConversionStatus.CONVERTED

            # check cxg metadata is updated
            array = tiledb.open(f"{testing_cxg_temp_directory}/cxg_group_metadata")
            actual_stored_metadata = dict(array.meta.items())

            expected_metadata = {
                "cxg_version": "1.0.0",
                "corpora": json.dumps(
                    {
                        "schema_version": "3.0.0",
                        "citation": "Publication www.doi.org/567.8",
                        "title": "New Dataset Title",
                    }
                ),
            }

            assert expected_metadata == actual_stored_metadata

    def test_update_cxg__with_spatial_deepzoom_assets(self):
        with tempfile.TemporaryDirectory() as tempdir:
            testing_cxg_temp_directory = tempdir

            cxg_metadata = {
                "cxg_version": "1.0.0",
                "corpora": json.dumps(
                    {"schema_version": "3.0.0", "citation": "Publication www.doi.org/123.4", "title": "Dataset Title 1"}
                ),
            }
            convert_dictionary_to_cxg_group(
                testing_cxg_temp_directory, cxg_metadata, group_metadata_name="cxg_group_metadata"
            )

            collection_version = self.generate_unpublished_collection(add_datasets=1)
            current_dataset_version = collection_version.datasets[0]
            new_dataset_version_id, _ = self.business_logic.ingest_dataset(
                collection_version_id=collection_version.version_id,
                url=None,
                file_size=0,
                current_dataset_version_id=current_dataset_version.version_id,
                start_step_function=False,
            )

            self.updater.s3_provider.upload_directory(
                tempdir,
                f"s3://{self.updater.spatial_deep_zoom_dir}/{current_dataset_version.version_id.id}/spatial.dzi",
            )

            metadata_update = DatasetArtifactMetadataUpdate(
                citation="Publication www.doi.org/567.8", title="New Dataset Title"
            )
            self.updater.update_cxg(
                None,
                testing_cxg_temp_directory,
                current_dataset_version.version_id,
                new_dataset_version_id,
                metadata_update,
            )

            # check new spatial deepzoom directory exists
            assert self.updater.s3_provider.uri_exists(
                f"s3://{self.updater.spatial_deep_zoom_dir}/{new_dataset_version_id.id}"
            )


class TestValidArtifactStatuses(BaseProcessingTest):
    def setUp(self):
        super().setUp()
        self.updater = DatasetMetadataUpdater(
            self.business_logic, "artifact_bucket", "cellxgene_bucket", "datasets_bucket", "spatial_deep_zoom_dir"
        )

    @parameterized.expand([DatasetConversionStatus.CONVERTED, DatasetConversionStatus.SKIPPED])
    def test_has_valid_artifact_statuses(self, rds_status):
        dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate(status_key=DatasetStatusKey.H5AD, status=DatasetConversionStatus.CONVERTED),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=rds_status),
                DatasetStatusUpdate(status_key=DatasetStatusKey.CXG, status=DatasetConversionStatus.CONVERTED),
            ]
        )

        dataset_version_id = DatasetVersionId(dataset_version.dataset_version_id)
        assert self.updater.has_valid_artifact_statuses(dataset_version_id) is True

    @parameterized.expand([DatasetConversionStatus.CONVERTING, DatasetConversionStatus.FAILED])
    def test_has_valid_artifact_statuses__invalid_rds_status(self, rds_status):
        dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate(status_key=DatasetStatusKey.H5AD, status=DatasetConversionStatus.CONVERTED),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=rds_status),
                DatasetStatusUpdate(status_key=DatasetStatusKey.CXG, status=DatasetConversionStatus.CONVERTED),
            ]
        )

        dataset_version_id = DatasetVersionId(dataset_version.dataset_version_id)
        assert self.updater.has_valid_artifact_statuses(dataset_version_id) is False

    @parameterized.expand([DatasetConversionStatus.CONVERTING, DatasetConversionStatus.FAILED])
    def test_has_valid_artifact_statuses__invalid_h5ad_status(self, h5ad_status):
        dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate(status_key=DatasetStatusKey.H5AD, status=h5ad_status),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
                DatasetStatusUpdate(status_key=DatasetStatusKey.CXG, status=DatasetConversionStatus.CONVERTED),
            ]
        )

        dataset_version_id = DatasetVersionId(dataset_version.dataset_version_id)
        assert self.updater.has_valid_artifact_statuses(dataset_version_id) is False

    @parameterized.expand([DatasetConversionStatus.CONVERTING, DatasetConversionStatus.FAILED])
    def test_has_valid_artifact_statuses__invalid_cxg_status(self, cxg_status):
        dataset_version = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate(status_key=DatasetStatusKey.H5AD, status=DatasetConversionStatus.CONVERTED),
                DatasetStatusUpdate(status_key=DatasetStatusKey.RDS, status=DatasetConversionStatus.CONVERTED),
                DatasetStatusUpdate(status_key=DatasetStatusKey.CXG, status=cxg_status),
            ]
        )

        dataset_version_id = DatasetVersionId(dataset_version.dataset_version_id)
        assert self.updater.has_valid_artifact_statuses(dataset_version_id) is False
