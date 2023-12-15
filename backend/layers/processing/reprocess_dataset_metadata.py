"""
Updates citations in-place across dataset artifacts for a Collection
"""
import json
import logging
import os
from multiprocessing import Process

import scanpy
import tiledb
from rpy2.robjects import ListVector, StrVector, r
from rpy2.robjects.packages import importr

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.exceptions import ProcessingFailed
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_download import ProcessDownload
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

base = importr("base")
seurat = importr("SeuratObject")

configure_logging(level=logging.INFO)


class DatasetMetadataReprocessWorker(ProcessDownload):
    def __init__(self, artifact_bucket: str, datasets_bucket: str) -> None:
        # init each worker with business logic backed by non-shared DB connection
        self.business_logic = BusinessLogic(
            DatabaseProvider(),
            None,
            None,
            None,
            S3Provider(),
            UriProvider(),
        )
        super().__init__(self.business_logic, self.business_logic.uri_provider, self.business_logic.s3_provider)
        self.artifact_bucket = artifact_bucket
        self.datasets_bucket = datasets_bucket

    def create_updated_artifacts(
        self,
        file_name: str,
        artifact_type: str,
        key_prefix: str,
        dataset_version_id: DatasetVersionId,
    ):
        try:
            s3_uri = self.upload_artifact(file_name, key_prefix, self.artifact_bucket)
            self.logger.info(f"Uploaded [{dataset_version_id}/{file_name}] to {s3_uri}")

            # TODO: include datasets_bucket uploads or not?
            # key = ".".join((key_prefix, DatasetArtifactType.H5AD))
            # self.s3_provider.upload_file(
            #     file_name, self.datasets_bucket, key, extra_args={"ACL": "bucket-owner-full-control"}
            # )
            # datasets_s3_uri = self.make_s3_uri(self.datasets_bucket, key_prefix, key)
            # self.logger.info(f"Uploaded [{dataset_version_id}/{file_name}] to {datasets_s3_uri}")
        except Exception:
            self.logger.error(f"Uploading Artifact {artifact_type} from dataset {dataset_version_id} failed.")
            raise ProcessingFailed from None

    def update_h5ad(
        self,
        h5ad_uri: str,
        current_dataset_version: DatasetVersion,
        new_key_prefix: str,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        h5ad_filename = self.download_from_source_uri(
            source_uri=h5ad_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )

        adata = scanpy.read_h5ad(h5ad_filename)
        metadata = current_dataset_version.metadata
        # maps artifact name for metadata field to DB field name, if different
        for key, val in metadata_update.as_dict_without_none_values().items():
            adata.uns[key] = val
            setattr(metadata, key, val)

        adata.write(h5ad_filename, compression="gzip")
        current_dataset_version_id = current_dataset_version.version_id
        self.business_logic.set_dataset_metadata(current_dataset_version_id, metadata)
        self.create_updated_artifacts(
            h5ad_filename, DatasetArtifactType.H5AD, new_key_prefix, current_dataset_version_id
        )
        os.remove(h5ad_filename)

    def update_rds(
        self,
        rds_uri: str,
        new_key_prefix: str,
        current_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        seurat_filename = self.download_from_source_uri(
            source_uri=rds_uri,
            local_path=CorporaConstants.LABELED_RDS_ARTIFACT_FILENAME,
        )

        rds_object = base.readRDS(seurat_filename)

        new_keys = []
        seurat_metadata = seurat.Misc(object=rds_object)
        for key, val in metadata_update.as_dict_without_none_values().items():
            if seurat_metadata.rx2[key]:
                val = val if isinstance(val, list) else [val]
                seurat_metadata[seurat_metadata.names.index(key)] = StrVector(val)
            else:
                new_keys.append((key, val))

        if new_keys:
            new_key_vector = ListVector({k: v for k, v in new_keys})
            seurat_metadata += new_key_vector
            r.assign("rds_object", rds_object)
            r.assign("seurat_metadata", seurat_metadata)
            r("rds_object@misc <- seurat_metadata")
            rds_object = r["rds_object"]

        base.saveRDS(rds_object, file=seurat_filename)

        self.create_updated_artifacts(
            seurat_filename, DatasetArtifactType.RDS, new_key_prefix, current_dataset_version_id
        )
        os.remove(seurat_filename)

    def update_cxg(
        self,
        cxg_uri: str,
        new_cxg_dir: str,
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


class DatasetMetadataReprocess(ProcessDownload):
    def __init__(
        self, business_logic: BusinessLogic, artifact_bucket: str, cellxgene_bucket: str, datasets_bucket: str
    ) -> None:
        super().__init__(business_logic, business_logic.uri_provider, business_logic.s3_provider)
        self.artifact_bucket = artifact_bucket
        self.cellxgene_bucket = cellxgene_bucket
        self.datasets_bucket = datasets_bucket

    @staticmethod
    def update_h5ad(
        artifact_bucket: str,
        datasets_bucket: str,
        h5ad_uri: str,
        current_dataset_version: DatasetVersion,
        new_key_prefix: str,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        DatasetMetadataReprocessWorker(artifact_bucket, datasets_bucket).update_h5ad(
            h5ad_uri,
            current_dataset_version,
            new_key_prefix,
            metadata_update,
        )

    @staticmethod
    def update_rds(
        artifact_bucket: str,
        datasets_bucket: str,
        rds_uri: str,
        new_key_prefix: str,
        current_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        DatasetMetadataReprocessWorker(artifact_bucket, datasets_bucket).update_rds(
            rds_uri, new_key_prefix, current_dataset_version_id, metadata_update
        )

    @staticmethod
    def update_cxg(
        artifact_bucket: str,
        datasets_bucket: str,
        cxg_uri: str,
        new_cxg_dir: str,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        DatasetMetadataReprocessWorker(artifact_bucket, datasets_bucket).update_cxg(
            cxg_uri, new_cxg_dir, metadata_update
        )

    def update_dataset_metadata(
        self,
        current_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        current_dataset_version = self.business_logic.get_dataset_version(current_dataset_version_id)
        if current_dataset_version.status.processing_status != DatasetProcessingStatus.SUCCESS:
            self.logger.info(
                f"Dataset {current_dataset_version_id} is not successfully processed. Skipping metadata update."
            )
            return

        artifact_uris = {artifact.type: artifact.uri for artifact in current_dataset_version.artifacts}

        new_artifact_key_prefix = self.get_key_prefix(current_dataset_version_id.id) + "_updated"

        artifact_jobs = []

        if DatasetArtifactType.H5AD in artifact_uris:
            self.logger.info("Main: Starting thread for h5ad update")
            h5ad_job = Process(
                target=DatasetMetadataReprocess.update_h5ad,
                args=(
                    self.artifact_bucket,
                    self.datasets_bucket,
                    artifact_uris[DatasetArtifactType.H5AD],
                    current_dataset_version,
                    new_artifact_key_prefix,
                    metadata_update,
                ),
            )
            artifact_jobs.append(h5ad_job)
            h5ad_job.start()
        else:
            self.logger.error(f"Cannot find labeled H5AD artifact uri for {current_dataset_version_id}.")
            raise ProcessingFailed from None

        if DatasetArtifactType.RDS in artifact_uris:
            self.logger.info("Main: Starting thread for rds update")
            # RDS-only, one-time
            metadata_update.schema_reference = (
                "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md"
            )
            rds_job = Process(
                target=DatasetMetadataReprocess.update_rds,
                args=(
                    self.artifact_bucket,
                    self.datasets_bucket,
                    artifact_uris[DatasetArtifactType.RDS],
                    new_artifact_key_prefix,
                    current_dataset_version_id,
                    metadata_update,
                ),
            )
            artifact_jobs.append(rds_job)
            rds_job.start()
        elif current_dataset_version.status.rds_status == DatasetConversionStatus.SKIPPED:
            pass
        else:
            self.logger.error(
                f"Cannot find RDS artifact uri for {current_dataset_version_id}, and Conversion Status is not SKIPPED."
            )
            raise ProcessingFailed from None

        if DatasetArtifactType.CXG in artifact_uris:
            self.logger.info("Main: Starting thread for cxg update")
            cxg_job = Process(
                target=DatasetMetadataReprocess.update_cxg,
                args=(
                    self.artifact_bucket,
                    self.datasets_bucket,
                    artifact_uris[DatasetArtifactType.CXG],
                    f"s3://{self.cellxgene_bucket}/{new_artifact_key_prefix}.cxg",
                    metadata_update,
                ),
            )
            artifact_jobs.append(cxg_job)
            cxg_job.start()
        else:
            self.logger.error(f"Cannot find cxg artifact uri for {current_dataset_version_id}.")
            raise ProcessingFailed from None

        # blocking call on async functions before checking for valid artifact statuses
        [j.join() for j in artifact_jobs]

    def update_dataset_citation(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
    ):

        collection_version = self.business_logic.get_collection_version(collection_version_id)
        doi = next((link.uri for link in collection_version.metadata.links if link.type == "DOI"), None)
        new_citation = self.business_logic.generate_dataset_citation(
            collection_version.collection_id, dataset_version_id, doi
        )
        metadata_update = DatasetArtifactMetadataUpdate(citation=new_citation)
        self.update_dataset_metadata(dataset_version_id, metadata_update)


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
    collection_version_id = CollectionVersionId(os.environ["COLLECTION_VERSION_ID"])
    dataset_version_id = DatasetVersionId(os.environ["DATASET_VERSION_ID"])
    DatasetMetadataReprocess(
        business_logic, artifact_bucket, cellxgene_bucket, datasets_bucket
    ).update_dataset_citation(collection_version_id, dataset_version_id)
