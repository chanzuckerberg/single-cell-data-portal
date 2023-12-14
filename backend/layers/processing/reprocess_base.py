"""
Base classes for reprocessing dataset artifacts.
"""
import logging
import os

import scanpy
from rpy2.robjects.packages import importr

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.process_download import ProcessDownload
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider

base = importr("base")
seurat = importr("SeuratObject")


class DatasetMetadataUpdaterWorker(ProcessDownload):
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
        self.logger = logging.getLogger(self.__class__.__name__)


class CXGWorkerBase(DatasetMetadataUpdaterWorker):
    def __init__(self, artifact_bucket: str, datasets_bucket: str, cxg_uri: str, new_cxg_dir: str) -> None:
        super().__init__(artifact_bucket, datasets_bucket)
        self.cxg_uri = cxg_uri
        self.new_cxg_dir = new_cxg_dir

    def run(self) -> None:
        self.s3_provider.upload_directory(self.cxg_uri, self.new_cxg_dir)
        self.update()
        self.post()

    def post(self) -> None:
        """
        Update the backend to reflect the changes
        :return:
        """
        pass

    def update(self) -> None:
        """
        Update the cxg object in the destination directory
        """
        pass


class H5ADWorkerBase(DatasetMetadataUpdaterWorker):
    def __init__(self, artifact_bucket: str, datasets_bucket: str, h5ad_uri: str) -> None:
        super().__init__(artifact_bucket, datasets_bucket)
        self.h5ad_uri = h5ad_uri

    def run(self) -> None:
        """
        Update the anndata object
        """
        h5ad_filename = self.download_from_source_uri(
            source_uri=self.h5ad_uri,
            local_path=CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME,
        )
        try:
            adata = scanpy.read_h5ad(h5ad_filename)
            self.update(adata)
            adata.write_h5ad(compression="gzip")
            self.post(h5ad_filename)
        finally:
            self.logger.info("Cleanup h5ad")
            os.remove(h5ad_filename)

    def post(self, h5ad_filename: str) -> None:
        """
        Update the backend to reflect the changes
        :param h5ad_filename:
        """
        pass

    def update(self, adata: scanpy.AnnData) -> None:
        """
        Update the anndata object in place
        :param adata:
        """
        pass


class RDSWorkerBase(DatasetMetadataUpdaterWorker):
    def __init__(self, artifact_bucket: str, datasets_bucket: str, rds_uri: str) -> None:
        super().__init__(artifact_bucket, datasets_bucket)
        self.rds_uri = rds_uri

    def run(self) -> None:
        """
        Update the seurat object
        """
        seurat_filename = self.download_from_source_uri(
            source_uri=self.rds_uri,
            local_path=CorporaConstants.LABELED_RDS_ARTIFACT_FILENAME,
        )
        try:
            rds_object = base.readRDS(seurat_filename)
            self.update(rds_object)
            base.saveRDS(rds_object, file=seurat_filename)
            self.post(seurat_filename)
        finally:
            self.logger.info("Cleanup seurat")
            os.remove(seurat_filename)

    def post(self, rds_filename: str) -> None:
        """
        Update the backend to reflect the changes
        :param rds_filename:
        """
        pass

    def update(self, rds_object) -> None:
        """
        Update the seurat object in place
        :param rds_object:
        """
        pass
