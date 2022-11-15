from backend.layers.processing.exceptions import ConversionFailed, ProcessingCanceled, ProcessingFailed, UploadFailed, ValidationFailed
from backend.layers.processing.process_logic import ProcessingLogic
import os
from backend.layers.processing.process_cxg import ProcessCxg

from backend.layers.processing.process_download_validate import ProcessDownloadValidate
from backend.layers.processing.process_seurat import ProcessSeurat
from entities import DatasetConversionStatus, DatasetProcessingStatus, DatasetStatusKey, DatasetUploadStatus, DatasetValidationStatus, DatasetVersionId


class ProcessMain(ProcessingLogic):

    process_download_validate: ProcessDownloadValidate
    process_seurat: ProcessSeurat
    process_cxg: ProcessCxg

    def log_batch_environment(self):
        batch_environment_variables = [
            "AWS_BATCH_CE_NAME",
            "AWS_BATCH_JOB_ATTEMPT",
            "AWS_BATCH_JOB_ID",
            "STEP_NAME",
            "DROPBOX_URL",
            "ARTIFACT_BUCKET",
            "CELLXGENE_BUCKET",
            "DATASET_ID",
            "DEPLOYMENT_STAGE",
            "MAX_ATTEMPTS",
        ]
        env_vars = dict()
        for var in batch_environment_variables:
            env_vars[var] = os.getenv(var)
        self.logger.info(f"Batch Job Info: {env_vars}")


    def process(self, dataset_id: DatasetVersionId, step_name: str, dropbox_uri: str, artifact_bucket: str, cxg_bucket: str):
        self.log_batch_environment()
        self.logger.info(f"Processing dataset {dataset_id}")
        try:
            if step_name == "download-validate":
                self.process_download_validate.process(dataset_id, dropbox_uri, artifact_bucket)
            elif step_name == "cxg":
                self.process_cxg.process(dataset_id, artifact_bucket, cxg_bucket)
            elif step_name == "seurat":
                self.process_seurat.process(dataset_id, artifact_bucket)
            elif step_name == "cxg_remaster":
                raise NotImplementedError("cxg remasters are not supported anymore")
            else:
                self.logger.error(f"Step function configuration error: Unexpected STEP_NAME '{step_name}'")

        # TODO: this could be better - maybe collapse all these exceptions and pass in the status key and value
        except ProcessingCanceled:
            pass # TODO: what's the effect of canceling a dataset now?
        except ValidationFailed:
            self.update_processing_status(dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID)
            return False
        except ProcessingFailed:
            self.update_processing_status(dataset_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE)
            return False
        except UploadFailed:
            self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.FAILED)
            return False
        except ConversionFailed as e:
            self.update_processing_status(dataset_id, e.failed_status, DatasetConversionStatus.FAILED)
            return False
        except Exception as e:
            self.logger.exception(f"An unexpected error occurred while processing the data set: {e}")
            if step_name == "download-validate":
                self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.FAILED)
            elif step_name == "seurat":
                self.update_processing_status(dataset_id, DatasetStatusKey.RDS, DatasetConversionStatus.FAILED)
            elif step_name == "cxg":
                self.update_processing_status(dataset_id, DatasetStatusKey.CXG, DatasetConversionStatus.FAILED)
            return False

        return True

    def main(self):
        dataset_id = os.environ["DATASET_ID"]
        step_name = os.environ["STEP_NAME"]
        dropbox_uri = os.environ["DROPBOX_URL"]
        artifact_bucket = os.environ["ARTIFACT_BUCKET"]
        cxg_bucket = os.environ["CELLXGENE_BUCKET"]
        rv = self.process(
            dataset_id=DatasetVersionId(dataset_id),
            step_name=step_name,
            dropbox_uri=dropbox_uri,
            artifact_bucket=artifact_bucket,
            cxg_bucket=cxg_bucket
        )
        return 0 if rv else 1

if __name__ == "__main__":
    pass # TODO: create all of the above