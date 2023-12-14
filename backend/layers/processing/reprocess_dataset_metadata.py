"""
Updates citations in-place across dataset artifacts for a Collection
"""
import json
import logging
import os
from multiprocessing import Process

import scanpy
import tiledb
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactMetadataUpdate,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.exceptions import ProcessingFailed
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_download import ProcessDownload
from backend.layers.processing.reprocess_base import CXGWorkerBase, H5ADWorkerBase, RDSWorkerBase
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

base = importr("base")
seurat = importr("SeuratObject")

configure_logging(level=logging.INFO)


class CreateUpdatedArtifactMixin:
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


class H5ADWorker(H5ADWorkerBase, CreateUpdatedArtifactMixin):
    def __init__(
        self,
        artifact_bucket: str,
        h5ad_uri: str,
        datasets_bucket: str,
        current_dataset_version: DatasetVersion,
        new_key_prefix: str,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        super().__init__(artifact_bucket, datasets_bucket, h5ad_uri)
        self.current_dataset_version = current_dataset_version
        self.new_key_prefix = new_key_prefix
        self.current_dataset_version = current_dataset_version
        self.metadata_update = metadata_update

    def act(self, adata: scanpy.AnnData):
        metadata = self.current_dataset_version.metadata
        # maps artifact name for metadata field to DB field name, if different
        for key, val in self.metadata_update.as_dict_without_none_values().items():
            adata.uns[key] = val
            setattr(metadata, key, val)

    def post(self, h5ad_filename: str):
        current_dataset_version_id = self.current_dataset_version.version_id
        self.business_logic.set_dataset_metadata(current_dataset_version_id, self.current_dataset_version.metadata)
        self.create_updated_artifacts(
            h5ad_filename, DatasetArtifactType.H5AD, self.new_key_prefix, current_dataset_version_id
        )


class RDSWorker(RDSWorkerBase, CreateUpdatedArtifactMixin):
    def __init__(
        self,
        artifact_bucket: str,
        datasets_bucket: str,
        rds_uri: str,
        new_key_prefix: str,
        current_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        super().__init__(artifact_bucket, datasets_bucket, rds_uri)
        self.new_key_prefix = new_key_prefix
        self.current_dataset_version_id = current_dataset_version_id
        self.metadata_update = metadata_update

    def act(self, rds_object):
        for key, val in self.metadata_update.as_dict_without_none_values().items():
            seurat_metadata = seurat.Misc(object=rds_object)
            if seurat_metadata.rx2[key]:
                val = val if isinstance(val, list) else [val]
                seurat_metadata[seurat_metadata.names.index(key)] = StrVector(val)

    def post(self, rds_filename: str):
        self.create_updated_artifacts(
            rds_filename, DatasetArtifactType.RDS, self.new_key_prefix, self.current_dataset_version_id
        )


class CXGWorker(CXGWorkerBase):
    def __init__(
        self,
        artifact_bucket: str,
        datasets_bucket: str,
        cxg_uri: str,
        new_cxg_dir: str,
        dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ) -> None:
        super().__init__(artifact_bucket, datasets_bucket, cxg_uri, new_cxg_dir)
        self.dataset_version_id = dataset_version_id
        self.metadata_update = metadata_update

    def act(self):
        ctx = tiledb.Ctx(H5ADDataFile.tile_db_ctx_config)
        array_name = f"{self.new_cxg_dir}/cxg_group_metadata"
        with tiledb.open(array_name, mode="r", ctx=ctx) as metadata_array:
            cxg_metadata_dict = json.loads(metadata_array.meta["corpora"])
            cxg_metadata_dict.update(self.metadata_update.as_dict_without_none_values())

        with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
            metadata_array.meta["corpora"] = json.dumps(cxg_metadata_dict)

    def post(self):
        self.business_logic.add_dataset_artifact(self.dataset_version_id, DatasetArtifactType.CXG, self.new_cxg_dir)
        self.update_processing_status(self.dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.CONVERTED)


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
        H5ADWorker(
            artifact_bucket, datasets_bucket, h5ad_uri, current_dataset_version, new_key_prefix, metadata_update
        ).update()

    @staticmethod
    def update_rds(
        artifact_bucket: str,
        datasets_bucket: str,
        rds_uri: str,
        new_key_prefix: str,
        current_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        RDSWorker(
            artifact_bucket, datasets_bucket, rds_uri, new_key_prefix, current_dataset_version_id, metadata_update
        ).update()

    @staticmethod
    def update_cxg(
        artifact_bucket: str,
        datasets_bucket: str,
        cxg_uri: str,
        new_cxg_dir: str,
        metadata_update: DatasetArtifactMetadataUpdate,
    ):
        CXGWorker(artifact_bucket, datasets_bucket, cxg_uri, new_cxg_dir, metadata_update).update()

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
