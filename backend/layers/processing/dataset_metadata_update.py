"""
Creates a new DatasetVersion to update metadata across dataset artifacts
"""
import json
import logging
import os
from typing import Dict

import scanpy

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_download import ProcessDownload
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

configure_logging(level=logging.INFO)


class DatasetMetadataUpdate(ProcessDownload):
    def __init__(
        self, business_logic: BusinessLogic, artifact_bucket: str, cellxgene_bucket: str, datasets_bucket: str
    ) -> None:
        super().__init__(business_logic, business_logic.uri_provider, business_logic.s3_provider)
        self.artifact_bucket = artifact_bucket
        self.cellxgene_bucket = cellxgene_bucket
        self.datasets_bucket = datasets_bucket

    def update_metadata(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ):
        # TODO: add error checks, skip datasets with non-success processing statuses
        original_dataset_version = self.business_logic.get_dataset_version(dataset_version_id)
        raw_h5ad_s3_uri = None
        h5ad_s3_uri = None
        rds_s3_uri = None
        cxg_s3_uri = None
        for artifact in original_dataset_version.artifacts:
            if artifact.type == DatasetArtifactType.RAW_H5AD:
                raw_h5ad_s3_uri = artifact.uri
            elif artifact.type == DatasetArtifactType.H5AD:
                h5ad_s3_uri = artifact.uri
            elif artifact.type == DatasetArtifactType.RDS:
                rds_s3_uri = artifact.uri
            elif artifact.type == DatasetArtifactType.CXG:
                cxg_s3_uri = artifact.uri
        # TODO: factor out the logic we want from ingest dataset from the URI check, so we don't need to pass a uri
        new_dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection_version_id, raw_h5ad_s3_uri, None, dataset_version_id, start_step_function=False
        )
        self.process(new_dataset_version_id, raw_h5ad_s3_uri, self.artifact_bucket)

        key_prefix = self.get_key_prefix(new_dataset_version_id)
        # update anndata
        # TODO: factor out any shared logic
        if h5ad_s3_uri:
            h5ad_filename = self.download_from_source_uri(
                source_uri=h5ad_s3_uri,
                local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
            )
            adata = scanpy.read_h5ad(h5ad_filename, backed="r")
            metadata = original_dataset_version.metadata
            for key, val in metadata_update_dict:
                adata.uns[key] = val
                if hasattr(metadata, key):
                    setattr(metadata, key, val)
            adata.write(h5ad_filename)
            self.business_logic.set_dataset_metadata(new_dataset_version_id, metadata)

            self.create_artifact(
                h5ad_filename,
                DatasetArtifactType.H5AD,
                key_prefix,
                new_dataset_version_id,
                self.artifact_bucket,
                DatasetStatusKey.H5AD,
                datasets_bucket=self.datasets_bucket,
            )
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED
            )
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID
            )

        # TODO: update seurat metadata dict
        if rds_s3_uri:
            seurat_filename = self.download_from_source_uri(
                source_uri=rds_s3_uri,
                local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,  # TODO: diff name?
            )
            self.create_artifact(
                seurat_filename,
                DatasetArtifactType.RDS,
                key_prefix,
                new_dataset_version_id,
                self.artifact_bucket,
                DatasetStatusKey.RDS,
                datasets_bucket=self.datasets_bucket,
            )
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.CONVERTED
            )

        if cxg_s3_uri:
            new_dir = f"s3://{cellxgene_bucket}/{key_prefix}.cxg/"
            self.s3_provider.upload_directory(cxg_s3_uri, new_dir)
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.CONVERTED
            )

        self.update_processing_status(
            new_dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
        )


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,  # Not required - decide if we should pass for safety
        None,  # Not required - decide if we should pass for safety
        S3Provider(),
        UriProvider(),
    )

    artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
    cellxgene_bucket = os.environ.get("CELLXGENE_BUCKET", "test-cellxgene-bucket")
    datasets_bucket = os.eniron.get("DATASETS_BUCKET", "test-datasets-bucket")
    collection_version_id = CollectionVersionId(os.environ["DATASET_VERSION_ID"])
    dataset_version_id = DatasetVersionId(os.environ["DATASET_VERSION_ID"])
    metadata_update_dict = json.loads(os.environ["METADATA_UPDATE_JSON"].strip())

    DatasetMetadataUpdate(business_logic, artifact_bucket, datasets_bucket).update_metadata(
        collection_version_id, dataset_version_id, metadata_update_dict
    )
