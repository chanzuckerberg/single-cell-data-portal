import os
import sys
from typing import Optional

from backend.layers.business.business import BusinessLogic
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.downloader import Downloader
from backend.layers.processing.exceptions import (
    ConversionFailed,
    ProcessingCanceled,
    ProcessingFailed,
    UploadFailed,
    ValidationFailed,
)
from backend.layers.processing.logger import configure_logging
from backend.layers.processing.process_cxg import ProcessCxg
from backend.layers.processing.process_download_validate import ProcessDownloadValidate
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.processing.process_seurat import ProcessSeurat
from backend.layers.processing.schema_migration import SchemaMigrate
from backend.layers.thirdparty.s3_provider import S3Provider, S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import (
    SchemaValidatorProvider,
    SchemaValidatorProviderInterface,
)
from backend.layers.thirdparty.uri_provider import UriProvider, UriProviderInterface

configure_logging()


class ProcessMain(ProcessingLogic):
    """
    Main class for the dataset pipeline processing
    """

    process_download_validate: ProcessDownloadValidate
    process_seurat: ProcessSeurat
    process_cxg: ProcessCxg

    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
        downloader: Downloader,
        schema_validator: SchemaValidatorProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider
        self.downloader = downloader
        self.schema_validator = schema_validator
        self.process_download_validate = ProcessDownloadValidate(
            self.business_logic, self.uri_provider, self.s3_provider, self.downloader, self.schema_validator
        )
        self.process_seurat = ProcessSeurat(self.business_logic, self.uri_provider, self.s3_provider)
        self.process_cxg = ProcessCxg(self.business_logic, self.uri_provider, self.s3_provider)
        self.schema_migrate = SchemaMigrate(business_logic)  # type: ignore

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
            "MIGRATE",
            "REMOTE_DEV_PREFIX",
        ]
        env_vars = dict()
        for var in batch_environment_variables:
            env_vars[var] = os.getenv(var)
        self.logger.info(f"Batch Job Info: {env_vars}")

    def process(
        self,
        dataset_id: DatasetVersionId,
        step_name: str,
        dropbox_uri: Optional[str],
        artifact_bucket: Optional[str],
        datasets_bucket: Optional[str],
        cxg_bucket: Optional[str],
    ):
        """
        Gets called by the step function at every different step, as defined by `step_name`
        """
        self.logger.info(f"Processing dataset {dataset_id}")  # type: ignore
        try:
            if step_name == "download-validate":
                self.process_download_validate.process(dataset_id, dropbox_uri, artifact_bucket, datasets_bucket)  # type: ignore
            elif step_name == "cxg":
                self.process_cxg.process(dataset_id, artifact_bucket, cxg_bucket)  # type: ignore
            elif step_name == "cxg_remaster":
                self.process_cxg.process(dataset_id, artifact_bucket, cxg_bucket, is_reprocess=True)  # type: ignore
            elif step_name == "seurat":
                self.process_seurat.process(dataset_id, artifact_bucket, datasets_bucket)  # type: ignore
            else:
                self.logger.error(f"Step function configuration error: Unexpected STEP_NAME '{step_name}'")  # type: ignore

        # TODO: this could be better - maybe collapse all these exceptions and pass in the status key and value
        except ProcessingCanceled:
            pass  # TODO: what's the effect of canceling a dataset now?
        except ValidationFailed as e:
            self.update_processing_status(
                dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID, e.errors
            )
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
            self.logger.exception(f"An unexpected error occurred while processing the data set: {e}")  # type: ignore
            if step_name == "download-validate":
                self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.FAILED)
            elif step_name == "seurat":
                self.update_processing_status(dataset_id, DatasetStatusKey.RDS, DatasetConversionStatus.FAILED)
            elif step_name == "cxg" or step_name == "cxg_remaster":
                self.update_processing_status(dataset_id, DatasetStatusKey.CXG, DatasetConversionStatus.FAILED)
            return False

        return True

    def main(self):
        self.log_batch_environment()
        step_name = os.environ["STEP_NAME"]
        if os.environ.get("MIGRATE"):
            self.logger.info("Migrating schema")
            rv = self.schema_migrate.migrate(step_name)
        else:
            dataset_id = os.environ["DATASET_ID"]

            dropbox_uri = os.environ.get("DROPBOX_URL")
            artifact_bucket = os.environ.get("ARTIFACT_BUCKET")
            datasets_bucket = os.environ.get("DATASETS_BUCKET")
            cxg_bucket = os.environ.get("CELLXGENE_BUCKET")
            rv = self.process(
                dataset_id=DatasetVersionId(dataset_id),
                step_name=step_name,
                dropbox_uri=dropbox_uri,
                artifact_bucket=artifact_bucket,
                datasets_bucket=datasets_bucket,
                cxg_bucket=cxg_bucket,
            )
        return 0 if rv else 1


if __name__ == "__main__":
    database_provider = DatabaseProvider()
    s3_provider = S3Provider()
    uri_provider = UriProvider()

    business_logic = BusinessLogic(
        database_provider,
        None,  # type: ignore# Not required - decide if we should pass for safety
        None,  # type: ignore# Not required - decide if we should pass for safety
        s3_provider,
        uri_provider,
    )

    downloader = Downloader(business_logic)
    schema_validator = SchemaValidatorProvider()

    process_main = ProcessMain(
        business_logic=business_logic,
        uri_provider=uri_provider,
        s3_provider=s3_provider,
        downloader=downloader,
        schema_validator=schema_validator,
    )

    rv = process_main.main()
    sys.exit(rv)
