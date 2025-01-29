"""
Creates a new DatasetVersion to update metadata across dataset artifacts
"""

import json
import logging
import os

import fsspec
import h5py
import tiledb

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.exceptions import ProcessingFailed
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_validate import ProcessValidate
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

configure_logging(level=logging.DEBUG)

# maps artifact name for metadata field to DB field name, if different
ARTIFACT_TO_DB_FIELD = {"title": "name"}
FIELDS_IN_RAW_H5AD = ["title"]


class DatasetMetadataUpdater(ProcessValidate):
    def __init__(
        self,
        business_logic: BusinessLogic,
        artifact_bucket: str,
        cellxgene_bucket: str,
        datasets_bucket: str,
        spatial_deep_zoom_dir: str,
    ) -> None:
        super().__init__(business_logic, business_logic.uri_provider, business_logic.s3_provider, None)
        self.artifact_bucket = artifact_bucket
        self.cellxgene_bucket = cellxgene_bucket
        self.datasets_bucket = datasets_bucket
        self.spatial_deep_zoom_dir = spatial_deep_zoom_dir
        self.fs = fsspec.filesystem("s3")

    def persist_artifact(
        self,
        src_uri: str,
        new_key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        filename: str,
        new_bucket: str,
        artifact_type: DatasetArtifactType,
    ) -> str:
        src_bucket, src_key = self.s3_provider.parse_s3_uri(src_uri)
        new_key = f"{new_key_prefix}/{filename}"
        self.s3_provider.copy_file(src_key=src_key, src_bucket=src_bucket, dst_key=new_key, dst_bucket=new_bucket)
        s3_uri = self.make_s3_uri(new_bucket, new_key_prefix, filename)
        self.business_logic.add_dataset_artifact(new_dataset_version_id, artifact_type, s3_uri)
        return s3_uri

    def update_raw_h5ad(
        self,
        raw_h5ad_uri: str,
        new_key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)
        new_s3_uri = self.persist_artifact(
            raw_h5ad_uri,
            new_key_prefix,
            new_dataset_version_id,
            CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
            self.artifact_bucket,
            DatasetArtifactType.RAW_H5AD,
        )
        s3_path = new_s3_uri.split("s3://")[-1]
        s3file = self.fs.open(s3_path)

        with h5py.File(s3file, "rb+") as f:
            for key, val in metadata_update.as_dict_without_none_values().items():
                if key in f["uns"]:
                    del f["uns"][key]
                f["uns"].create_dataset(key, data=val)

        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)

    def update_h5ad(
        self,
        h5ad_uri: str,
        current_dataset_version: DatasetVersion,
        new_key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTING)
        new_s3_uri = self.persist_artifact(
            h5ad_uri,
            new_key_prefix,
            new_dataset_version_id,
            CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
            self.artifact_bucket,
            DatasetArtifactType.H5AD,
        )
        s3_path = new_s3_uri.split("s3://")[-1]
        s3file = self.fs.open(s3_path)

        metadata = current_dataset_version.metadata
        # maps artifact name for metadata field to DB field name, if different
        with h5py.File(s3file, "rb+") as f:
            for key, val in metadata_update.as_dict_without_none_values().items():
                if key in f["uns"]:
                    del f["uns"][key]
                f["uns"].create_dataset(key, data=val)

                db_field = ARTIFACT_TO_DB_FIELD.get(key) if key in ARTIFACT_TO_DB_FIELD else key
                setattr(metadata, db_field, val)

        self.business_logic.set_dataset_metadata(new_dataset_version_id, metadata)

        new_key = f"{new_key_prefix}/{CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME}"
        public_bucket_key = f"{new_dataset_version_id}.h5ad"
        self.s3_provider.copy_file(
            src_key=new_key,
            src_bucket=self.artifact_bucket,
            dst_key=public_bucket_key,
            dst_bucket=self.datasets_bucket,
        )

        self.update_processing_status(
            new_dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)

    def update_cxg(
        self,
        cxg_uri: str,
        new_cxg_dir: str,
        current_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        self.s3_provider.upload_directory(cxg_uri, new_cxg_dir)

        current_spatial_deep_zoom_dir = f"s3://{self.spatial_deep_zoom_dir}/{current_dataset_version_id.id}"
        spatial_dzi_uri = f"{current_spatial_deep_zoom_dir}/spatial.dzi"
        # Copy spatial deep zoom directory if it exists (only exists for Visium datasets)
        if self.s3_provider.uri_exists(spatial_dzi_uri):
            new_spatial_deep_zoom_dir = f"s3://{self.spatial_deep_zoom_dir}/{new_dataset_version_id.id}"
            self.s3_provider.upload_directory(current_spatial_deep_zoom_dir, new_spatial_deep_zoom_dir)

        ctx = tiledb.Ctx(H5ADDataFile.tile_db_ctx_config)
        array_name = f"{new_cxg_dir}/cxg_group_metadata"
        with tiledb.open(array_name, mode="r", ctx=ctx) as metadata_array:
            cxg_metadata_dict = json.loads(metadata_array.meta["corpora"])
            cxg_metadata_dict.update(metadata_update.as_dict_without_none_values())

        with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
            metadata_array.meta["corpora"] = json.dumps(cxg_metadata_dict)

        self.business_logic.add_dataset_artifact(new_dataset_version_id, DatasetArtifactType.CXG, new_cxg_dir)
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.CONVERTED)

    def update_metadata(
        self,
        current_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        current_dataset_version = self.business_logic.get_dataset_version(current_dataset_version_id)
        if current_dataset_version.status.processing_status != DatasetProcessingStatus.SUCCESS:
            self.logger.info(
                f"Dataset {current_dataset_version_id} is not successfully processed. Skipping metadata update."
            )
            return

        artifact_uris = {artifact.type: artifact.uri for artifact in current_dataset_version.artifacts}

        new_dataset_version = self.business_logic.get_dataset_version(new_dataset_version_id)
        if new_dataset_version.status.processing_status != DatasetProcessingStatus.INITIALIZED:
            self.logger.info(
                f"Dataset {new_dataset_version_id} has processing status "
                f"{new_dataset_version.status.processing_status} rather than expected INITIALIZED."
                f"Skipping metadata update."
            )
            return

        new_artifact_key_prefix = self.get_key_prefix(new_dataset_version_id.id)
        if DatasetArtifactType.RAW_H5AD in artifact_uris:
            raw_h5ad_uri = artifact_uris[DatasetArtifactType.RAW_H5AD]
        else:
            self.logger.error(f"Cannot find raw H5AD artifact uri for {current_dataset_version_id}.")
            raise ValueError

        # Only trigger raw H5AD update if any updated metadata is part of the raw H5AD artifact
        if any(getattr(metadata_update, field, None) for field in FIELDS_IN_RAW_H5AD):
            self.logger.info("Raw h5ad update required")
            self.update_raw_h5ad(
                raw_h5ad_uri,
                new_artifact_key_prefix,
                new_dataset_version_id,
                metadata_update,
            )
        else:
            self.logger.info("No raw h5ad update required")
            self.persist_artifact(
                raw_h5ad_uri,
                new_artifact_key_prefix,
                new_dataset_version_id,
                CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
                self.artifact_bucket,
                DatasetArtifactType.RAW_H5AD,
            )
            self.update_processing_status(new_dataset_version_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)

        if DatasetArtifactType.H5AD in artifact_uris:
            self.logger.info("Starting labeled h5ad update")
            self.update_h5ad(
                artifact_uris[DatasetArtifactType.H5AD],
                current_dataset_version,
                new_artifact_key_prefix,
                new_dataset_version_id,
                metadata_update,
            )
        else:
            self.logger.error(f"Cannot find labeled H5AD artifact uri for {current_dataset_version_id}.")
            self.update_processing_status(new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.FAILED)

        # Mark all RDS conversions as skipped
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED)

        if DatasetArtifactType.CXG in artifact_uris:
            self.logger.info("Starting cxg update")
            self.update_cxg(
                artifact_uris[DatasetArtifactType.CXG],
                f"s3://{self.cellxgene_bucket}/{new_artifact_key_prefix}.cxg",
                current_dataset_version_id,
                new_dataset_version_id,
                metadata_update,
            )
        else:
            self.logger.error(f"Cannot find cxg artifact uri for {current_dataset_version_id}.")
            self.update_processing_status(new_dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.FAILED)

        if self.has_valid_artifact_statuses(new_dataset_version_id):
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
            )
        else:
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE
            )
            status = self.business_logic.get_dataset_version(new_dataset_version_id).status
            raise ProcessingFailed(f"Artifact reprocessing failed, with statuses: {status.to_dict()}")

    def has_valid_artifact_statuses(self, dataset_version_id: DatasetVersionId) -> bool:
        dataset_version = self.business_logic.get_dataset_version(dataset_version_id)
        return (
            dataset_version.status.h5ad_status == DatasetConversionStatus.CONVERTED
            and dataset_version.status.cxg_status == DatasetConversionStatus.CONVERTED
            and (
                dataset_version.status.rds_status == DatasetConversionStatus.CONVERTED
                or dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED
            )
        )


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,
        None,
        None,
        S3Provider(),
        UriProvider(),
    )

    artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
    cellxgene_bucket = os.environ.get("CELLXGENE_BUCKET", "test-cellxgene-bucket")
    datasets_bucket = os.environ.get("DATASETS_BUCKET", "test-datasets-bucket")
    spatial_deep_zoom_bucket = os.environ.get("SPATIAL_DEEP_ZOOM_BUCKET", "test-spatial-deep-zoom-bucket")
    spatial_deep_zoom_dir = f"{spatial_deep_zoom_bucket}/spatial-deep-zoom"
    current_dataset_version_id = DatasetVersionId(os.environ["CURRENT_DATASET_VERSION_ID"])
    new_dataset_version_id = DatasetVersionId(os.environ["NEW_DATASET_VERSION_ID"])
    metadata_update = DatasetArtifactMetadataUpdate(**json.loads(os.environ["METADATA_UPDATE_JSON"]))
    DatasetMetadataUpdater(
        business_logic, artifact_bucket, cellxgene_bucket, datasets_bucket, spatial_deep_zoom_dir
    ).update_metadata(current_dataset_version_id, new_dataset_version_id, metadata_update)
