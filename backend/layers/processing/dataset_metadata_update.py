"""
Creates a new DatasetVersion to update metadata across dataset artifacts
"""
import json
import logging
import multiprocessing
import os

import scanpy
import tiledb
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_download import ProcessDownload
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

base = importr("base")
seurat = importr("SeuratObject")

configure_logging(level=logging.INFO)

# maps artifact name for metadata field to DB field name, if different
ARTIFACT_TO_DB_FIELD = {"title": "name"}


class DatasetMetadataUpdater(ProcessDownload):
    def __init__(
        self, business_logic: BusinessLogic, artifact_bucket: str, cellxgene_bucket: str, datasets_bucket: str
    ) -> None:
        super().__init__(business_logic, business_logic.uri_provider, business_logic.s3_provider)
        self.artifact_bucket = artifact_bucket
        self.cellxgene_bucket = cellxgene_bucket
        self.datasets_bucket = datasets_bucket

    def update_h5ad(
        self,
        h5ad_uri: str,
        current_dataset_version: DatasetVersion,
        new_key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        h5ad_filename = self.download_from_source_uri(
            source_uri=h5ad_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )

        adata = scanpy.read_h5ad(h5ad_filename, backed="r")
        metadata = current_dataset_version.metadata
        # maps artifact name for metadata field to DB field name, if different
        for key, val in metadata_update.as_dict_without_none_values().items():
            adata.uns[key] = val

            db_field = ARTIFACT_TO_DB_FIELD.get(key) if key in ARTIFACT_TO_DB_FIELD else key
            setattr(metadata, db_field, val)

        adata.write(h5ad_filename, compression="gzip")
        self.business_logic.set_dataset_metadata(new_dataset_version_id, metadata)

        self.create_artifact(
            h5ad_filename,
            DatasetArtifactType.H5AD,
            new_key_prefix,
            new_dataset_version_id,
            self.artifact_bucket,
            DatasetStatusKey.H5AD,
            datasets_bucket=self.datasets_bucket,
        )
        os.remove(h5ad_filename)
        self.update_processing_status(
            new_dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)

    def update_rds(
        self,
        rds_uri: str,
        new_key_prefix: str,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        seurat_filename = self.download_from_source_uri(
            source_uri=rds_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.CONVERTING)

        rds_object = base.readRDS(seurat_filename)

        for key, val in metadata_update.as_dict_without_none_values().items():
            seurat_metadata = seurat.Misc(object=rds_object)
            if seurat_metadata.rx2[key]:
                val = val if isinstance(val, list) else [val]
                seurat_metadata[seurat_metadata.names.index(key)] = StrVector(val)

        base.saveRDS(rds_object, file=seurat_filename)

        self.create_artifact(
            seurat_filename,
            DatasetArtifactType.RDS,
            new_key_prefix,
            new_dataset_version_id,
            self.artifact_bucket,
            DatasetStatusKey.RDS,
            datasets_bucket=self.datasets_bucket,
        )
        os.remove(seurat_filename)
        self.update_processing_status(new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.CONVERTED)

    def update_cxg(
        self,
        cxg_uri: str,
        new_cxg_dir: str,
        dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        self.s3_provider.upload_directory(cxg_uri, new_cxg_dir)
        ctx = tiledb.Ctx(H5ADDataFile.tile_db_ctx_config)
        array_name = f"{new_cxg_dir}/cxg_group_metadata"
        with tiledb.open(array_name, mode="r", ctx=ctx) as metadata_array:
            cxg_metadata_dict = json.loads(metadata_array.meta["corpora"])
            cxg_metadata_dict.update(metadata_update.as_dict_without_none_values())

        with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
            metadata_array.meta["corpora"] = json.dumps(cxg_metadata_dict)

        self.business_logic.add_dataset_artifact(dataset_version_id, DatasetArtifactType.CXG, new_cxg_dir)
        self.update_processing_status(dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.CONVERTED)

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

        if DatasetArtifactType.RAW_H5AD in artifact_uris:
            raw_h5ad_uri = artifact_uris[DatasetArtifactType.RAW_H5AD]
        else:
            self.logger.error(f"Cannot find raw H5AD artifact uri for {current_dataset_version_id}.")
            raise ValueError

        self.process(new_dataset_version_id, raw_h5ad_uri, self.artifact_bucket)

        new_artifact_key_prefix = self.get_key_prefix(new_dataset_version_id.id)

        artifact_jobs = []

        with multiprocessing.Pool() as pool:
            if DatasetArtifactType.H5AD in artifact_uris:
                self.logger.info("Main: Starting thread for h5ad update")
                artifact_jobs.append(
                    pool.apply_async(
                        self.update_h5ad,
                        (
                            artifact_uris[DatasetArtifactType.H5AD],
                            current_dataset_version,
                            new_artifact_key_prefix,
                            new_dataset_version_id,
                            metadata_update,
                        ),
                    )
                )
            else:
                self.logger.error(f"Cannot find labeled H5AD artifact uri for {current_dataset_version_id}.")
                self.update_processing_status(
                    new_dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.FAILED
                )

            if DatasetArtifactType.RDS in artifact_uris:
                self.logger.info("Main: Starting thread for rds update")
                artifact_jobs.append(
                    pool.apply_async(
                        self.update_rds,
                        (
                            artifact_uris[DatasetArtifactType.RDS],
                            new_artifact_key_prefix,
                            new_dataset_version_id,
                            metadata_update,
                        ),
                    )
                )
            elif current_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED:
                self.update_processing_status(
                    new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED
                )
            else:
                self.logger.error(
                    f"Cannot find RDS artifact uri for {current_dataset_version_id}, and Conversion Status is not SKIPPED."
                )
                self.update_processing_status(
                    new_dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.FAILED
                )

            if DatasetArtifactType.CXG in artifact_uris:
                self.logger.info("Main: Starting thread for cxg update")
                artifact_jobs.append(
                    pool.apply_async(
                        self.update_cxg,
                        (
                            artifact_uris[DatasetArtifactType.CXG],
                            f"s3://{self.cellxgene_bucket}/{new_artifact_key_prefix}.cxg",
                            new_dataset_version_id,
                            metadata_update,
                        ),
                    )
                )
            else:
                self.logger.error(f"Cannot find cxg artifact uri for {current_dataset_version_id}.")
                self.update_processing_status(
                    new_dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.FAILED
                )

        # blocking call on async jobs before checking for valid artifact statuses
        [j.get() for j in artifact_jobs]

        if self.has_valid_artifact_statuses(new_dataset_version_id):
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
            )
        else:
            self.update_processing_status(
                new_dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE
            )

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
        S3Provider(),
        UriProvider(),
    )

    artifact_bucket = os.environ.get("ARTIFACT_BUCKET", "test-bucket")
    cellxgene_bucket = os.environ.get("CELLXGENE_BUCKET", "test-cellxgene-bucket")
    datasets_bucket = os.environ.get("DATASETS_BUCKET", "test-datasets-bucket")
    current_dataset_version_id = DatasetVersionId(os.environ["CURRENT_DATASET_VERSION_ID"])
    new_dataset_version_id = DatasetVersionId(os.environ["NEW_DATASET_VERSION_ID"])
    metadata_update = DatasetArtifactMetadataUpdate(**json.loads(os.environ["METADATA_UPDATE_JSON"]))
    DatasetMetadataUpdater(business_logic, artifact_bucket, cellxgene_bucket, datasets_bucket).update_metadata(
        current_dataset_version_id, new_dataset_version_id, metadata_update
    )
