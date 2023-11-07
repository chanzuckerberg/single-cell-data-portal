"""
Creates a new DatasetVersion to update metadata across dataset artifacts
"""
import json
import logging
import os
from typing import Dict

import scanpy
import tiledb
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr

from backend.common.utils.corpora_constants import CorporaConstants
from backend.common.utils.tiledb import consolidation_buffer_size
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersion,
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

    def update_h5ad(
        self,
        h5ad_s3_uri: str,
        original_dataset_version: DatasetVersion,
        key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ):
        h5ad_filename = self.download_from_source_uri(
            source_uri=h5ad_s3_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )

        adata = scanpy.read_h5ad(h5ad_filename, backed="r")
        metadata = original_dataset_version.metadata
        for key, val in metadata_update_dict.items():
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
            new_dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)

    def update_rds(
        self,
        rds_s3_uri: str,
        key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ):
        seurat_filename = self.download_from_source_uri(
            source_uri=rds_s3_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.CONVERTING)

        base = importr("base")
        seurat = importr("SeuratObject")

        rds_object = base.readRDS(seurat_filename)

        for key, val in metadata_update_dict.items():
            seurat_metadata = seurat.Misc(object=rds_object)
            if seurat_metadata.rx2[key]:
                seurat_metadata[seurat_metadata.names.index(key)] = StrVector(list(val))

        base.saveRDS(rds_object, file=seurat_filename)

        self.create_artifact(
            seurat_filename,
            DatasetArtifactType.RDS,
            key_prefix,
            new_dataset_version_id,
            self.artifact_bucket,
            DatasetStatusKey.RDS,
            datasets_bucket=self.datasets_bucket,
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.CONVERTED)

    def update_cxg(
        self,
        cxg_s3_uri: str,
        key_prefix: str,
        dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ):
        new_cxg_dir = f"s3://{self.cellxgene_bucket}/{key_prefix}.cxg/"
        self.s3_provider.upload_directory(cxg_s3_uri, new_cxg_dir)
        ctx = tiledb.Ctx(
            {
                "sm.consolidation.buffer_size": consolidation_buffer_size(0.1),
                "py.deduplicate": True,  # May reduce memory requirements at cost of performance
            }
        )
        array_name = f"{new_cxg_dir}/cxg_group_metadata"
        with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
            for key, value in metadata_update_dict.items():
                metadata_array.meta[key] = value
        self.update_processing_status(dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.CONVERTED)

    def update_metadata(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ):
        original_dataset_version = self.business_logic.get_dataset_version(dataset_version_id)
        if original_dataset_version.status.processing_status != DatasetProcessingStatus.SUCCESS:
            self.logger.info(f"Dataset {dataset_version_id} is not successfully processed. Skipping metadata update.")
            return

        artifact_uris = {artifact.type: artifact.uri for artifact in original_dataset_version.artifacts}
        raw_h5ad_uri = artifact_uris[DatasetArtifactType.RAW_H5AD]

        new_dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection_version_id=collection_version_id,
            url=raw_h5ad_uri,
            file_size=0,
            existing_dataset_version_id=dataset_version_id,
            start_step_function=False,
        )

        self.process(new_dataset_version_id, raw_h5ad_uri, self.artifact_bucket)

        key_prefix = self.get_key_prefix(new_dataset_version_id.id)

        if DatasetArtifactType.H5AD in artifact_uris:
            self.update_h5ad(
                artifact_uris[DatasetArtifactType.H5AD],
                original_dataset_version,
                key_prefix,
                new_dataset_version_id,
                metadata_update_dict,
            )
        else:
            self.logger.error(f"Cannot find labeled H5AD artifact uri for {dataset_version_id}.")
            raise ValueError

        if DatasetArtifactType.RDS in artifact_uris:
            self.update_rds(
                artifact_uris[DatasetArtifactType.RDS], key_prefix, new_dataset_version_id, metadata_update_dict
            )
        elif original_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED:
            self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED)
        else:
            self.logger.error(
                f"Cannot find RDS artifact uri for {dataset_version_id}, " f"and Conversion Status is not SKIPPED."
            )
            raise ValueError

        if DatasetArtifactType.CXG in artifact_uris:
            self.update_cxg(
                artifact_uris[DatasetArtifactType.CXG], key_prefix, new_dataset_version_id, metadata_update_dict
            )
        else:
            self.logger.error(f"Cannot find cxg artifact uri for {dataset_version_id}.")
            raise ValueError

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
    datasets_bucket = os.environ.get("DATASETS_BUCKET", "test-datasets-bucket")
    collection_version_id = CollectionVersionId(os.environ["COLLECTION_VERSION_ID"])
    dataset_version_id = DatasetVersionId(os.environ["DATASET_VERSION_ID"])
    metadata_update_dict = json.loads(os.environ["METADATA_UPDATE_JSON"])

    DatasetMetadataUpdate(business_logic, artifact_bucket, cellxgene_bucket, datasets_bucket).update_metadata(
        collection_version_id, dataset_version_id, metadata_update_dict
    )
