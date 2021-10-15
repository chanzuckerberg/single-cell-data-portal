#!/usr/bin/env python3


"""
# Processing Status
## Initial
The initial processing_status when the container first runs is:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.WAITING
    upload_progress = 0
    upload_message = ""
    validation_status = ValidationStatus
    validation_message = ""
    rds_status = ConversionStatus
    cxg_status = ConversionStatus
    h5ad_status = ConversionStatus
}

## Upload

While uploading, upload_status changes UploadStatus.UPLOADING and upload_progress is updated regularly.
The processing_status should look like this:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADING
    upload_progress = 0.25
}

If upload succeeds the processing_status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
}


If upload fails the processing_status changes to:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.FAILED
    upload_progress = 0.25
    upload_message = "Some message"
}

## Validation
After upload, validation starts and processing status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALIDATING
}

If validation succeeds the process_status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALID
    rds_status = ConversionStatus.CONVERTING
    cxg_status = ConversionStatus.CONVERTING
    h5ad_status = ConversionStatus.CONVERTING
}

If validation fails the processing_status change to:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.FAILED,
    h5ad_status = ConversionStatus.CONVERTED
}

## Conversion
As each conversion is started the status for the file format is set to CONVERTING then CONVERTED if the conversion
completes successfully. The status is then set to UPLOADING while the file is copied to s3 and UPLOADED on success.
 Cellxgene data is converted/uploaded first. The h5ad_status is set to CONVERTED after validation/label writing is
 complete.
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    rds_status = ConversionStatus.CONVERTING
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.CONVERTED
}

If a conversion fails the processing_status will indicated it as follow:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    rds_status = ConversionStatus.FAILED
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.UPLOADED
}

Once all conversion are complete, the conversion status for each file will be either UPLOADED or FAILED:
{
    processing_status = ProcessingStatus.SUCCESS
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    rds_status = ConversionStatus.FAILED
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.UPLOADED
}
"""

import logging
import os
import subprocess
import sys
import typing
from datetime import datetime
from os.path import join

import numpy
import requests
import scanpy

from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import (
    ConversionStatus,
    DatasetArtifactFileType,
    DatasetArtifactType,
    ProcessingStatus,
    ValidationStatus,
)
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.dl_sources.url import from_url
from backend.corpora.dataset_processing.download import download
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ProcessingFailed, ValidationFailed
from backend.corpora.dataset_processing.h5ad_data_file import H5ADDataFile
from backend.corpora.dataset_processing.slack import format_slack_message

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# This is unfortunate, but this information doesn't appear to live anywhere
# accessible to the uploader
DEPLOYMENT_STAGE_TO_URL = {
    "staging": "https://cellxgene.staging.single-cell.czi.technology/e",
    "prod": "https://cellxgene.cziscience.com/e",
    "rdev": os.environ.get("FRONTEND_URL"),
    "dev": "https://cellxgene.dev.single-cell.czi.technology/e",
    "test": "http://frontend.corporanet.local:3000",
}


def check_env():
    """Verify that the required environment variables are set."""

    missing = []
    for env_var in ["DROPBOX_URL", "ARTIFACT_BUCKET", "CELLXGENE_BUCKET", "DATASET_ID", "DEPLOYMENT_STAGE"]:
        if env_var not in os.environ:
            missing.append(env_var)
    if missing:
        raise EnvironmentError(f"Missing environment variables: {missing}")


def create_artifact(
    file_name: str,
    artifact_type: DatasetArtifactFileType,
    bucket_prefix: str,
    dataset_id: str,
    artifact_bucket: str,
    processing_status_type: str,
):
    update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.UPLOADING})
    logger.info(f"Uploading [{dataset_id}/{file_name}] to S3 bucket: [{artifact_bucket}].")
    try:
        s3_uri = DatasetAsset.upload(file_name, bucket_prefix, artifact_bucket)
        with db_session_manager() as session:
            logger.info(f"Updating database with  {artifact_type}.")
            DatasetAsset.create(
                session,
                dataset_id=dataset_id,
                filename=file_name,
                filetype=artifact_type,
                type_enum=DatasetArtifactType.REMIX,
                user_submitted=True,
                s3_uri=s3_uri,
            )
        update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.UPLOADED})

    except Exception as e:
        update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.FAILED})
        raise e


def create_artifacts(local_filename, dataset_id, artifact_bucket):
    bucket_prefix = get_bucket_prefix(dataset_id)
    logger.info("Creating Artifacts.")
    # upload AnnData
    create_artifact(
        local_filename, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )

    # Process and upload seurat
    seurat_filename = convert_file_ignore_exceptions(
        make_seurat, local_filename, "Issue creating seurat.", dataset_id, "rds_status"
    )
    if seurat_filename:
        create_artifact(
            seurat_filename, DatasetArtifactFileType.RDS, bucket_prefix, dataset_id, artifact_bucket, "rds_status"
        )


def cancel_dataset(dataset_id):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        dataset.asset_deletion()
        dataset.delete()
        logger.info("Upload Canceled.")


def update_db(dataset_id, metadata=None, processing_status=None):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        if dataset.tombstone:
            raise ProcessingCancelled

        if metadata:
            logger.debug("Updating metadata.")
            dataset.update(**metadata)

        if processing_status:
            logger.debug(f"updating processing_status.{processing_status}")
            processing_status_updater(session, dataset.processing_status.id, processing_status)


def download_from_dropbox_url(dataset_uuid: str, dropbox_url: str, local_path: str) -> str:
    """Given a dropbox url, download it to local_path.
    Handles fixing the url so it downloads directly.
    """

    file_url = from_url(dropbox_url)
    if not file_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    file_info = file_url.file_info()
    status = download(dataset_uuid, file_url.url, local_path, file_info["size"])
    logger.info(status)
    logger.info("Download complete")
    return local_path


def extract_metadata(filename):
    """Pull metadata out of the AnnData file to insert into the dataset table."""

    adata = scanpy.read_h5ad(filename, backed="r")

    # TODO: Concern with respect to previous use of raising error when there is no raw layer.
    # This new way defaults to adata.X.
    if adata.raw is not None and adata.raw.X is not None:
        layer_for_mean_genes_per_cell = adata.raw.X
    else:
        layer_for_mean_genes_per_cell = adata.X

    # For mean_genes_per_cell, we only want the columns (genes) that have a feature_biotype of `gene`,
    # as opposed to `spike-in`
    filter_gene_vars = numpy.where(adata.var.feature_biotype == "gene")[0]

    # Calling np.count_nonzero on and h5py.Dataset appears to read the entire thing
    # into memory, so we need to chunk it to be safe.
    stride = 50000
    numerator, denominator = 0, 0
    for bounds in zip(
        range(0, layer_for_mean_genes_per_cell.shape[0], stride),
        range(stride, layer_for_mean_genes_per_cell.shape[0] + stride, stride),
    ):
        chunk = layer_for_mean_genes_per_cell[bounds[0] : bounds[1], filter_gene_vars]
        numerator += chunk.nnz if hasattr(chunk, "nnz") else numpy.count_nonzero(chunk)
        denominator += chunk.shape[0]

    def _get_term_pairs(base_term):
        base_term_id = base_term + "_ontology_term_id"
        return [
            {"label": k[0], "ontology_term_id": k[1]}
            for k in adata.obs.groupby([base_term, base_term_id]).groups.keys()
        ]

    def _get_is_primary_data():
        is_primary_data = adata.obs["is_primary_data"]
        if all(is_primary_data):
            return "PRIMARY"
        elif not any(is_primary_data):
            return "SECONDARY"
        else:
            return "BOTH"

    def _get_x_approximate_distribution():
        if "X_approximate_distribution" in adata.uns:
            return adata.uns["X_approximate_distribution"].upper()
        else:
            return None

    metadata = {
        "name": adata.uns["title"],
        "organism": _get_term_pairs("organism"),
        "tissue": _get_term_pairs("tissue"),
        "assay": _get_term_pairs("assay"),
        "disease": _get_term_pairs("disease"),
        "sex": _get_term_pairs("sex"),
        "ethnicity": _get_term_pairs("ethnicity"),
        "development_stage": _get_term_pairs("development_stage"),
        "cell_count": adata.shape[0],
        "mean_genes_per_cell": numerator / denominator,
        "is_primary_data": _get_is_primary_data(),
        "cell_type": _get_term_pairs("cell_type"),
        "x_normalization": adata.uns["X_normalization"],
        "x_approximate_distribution": _get_x_approximate_distribution(),
        "schema_version": adata.uns["schema_version"],
    }
    logger.info(f"Extract metadata: {metadata}")
    return metadata


def make_seurat(local_filename):
    """Create a Seurat rds file from the AnnData file."""

    try:
        subprocess.run(
            ["Rscript", os.path.join(os.path.abspath(os.path.dirname(__file__)), "make_seurat.R"), local_filename],
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError as ex:
        msg = f"Seurat conversion failed: {ex.output} {ex.stderr}"
        logger.exception(msg)
        raise RuntimeError(msg) from ex

    return local_filename.replace(".h5ad", ".rds")


def make_cxg(local_filename):
    """
    Convert the uploaded H5AD file to the CXG format servicing the cellxgene Explorer.
    """

    cxg_output_container = local_filename.replace(".h5ad", ".cxg")
    try:
        h5ad_data_file = H5ADDataFile(local_filename, var_index_column_name="feature_name")
        h5ad_data_file.to_cxg(cxg_output_container, 10.0)
    except Exception as ex:
        msg = "CXG conversion failed."
        logger.exception(msg)
        raise RuntimeError(msg) from ex

    return cxg_output_container


def copy_cxg_files_to_cxg_bucket(cxg_dir, object_key, cellxgene_bucket):
    """
    Copy cxg files to the cellxgene bucket (under the given object key) for access by the explorer
    returns the s3_uri where the cxg is stored
    """
    command = ["aws"]
    s3_uri = f"s3://{cellxgene_bucket}/{object_key}.cxg/"
    if os.getenv("BOTO_ENDPOINT_URL"):
        command.append(f"--endpoint-url={os.getenv('BOTO_ENDPOINT_URL')}")

    command.extend(
        [
            "s3",
            "cp",
            cxg_dir,
            s3_uri,
            "--recursive",
            "--acl",
            "bucket-owner-full-control",
        ]
    )
    subprocess.run(
        command,
        check=True,
    )
    return s3_uri


def convert_file_ignore_exceptions(
    converter: typing.Callable, local_filename: str, error_message: str, dataset_id: str, processing_status_type: str
) -> str:
    logger.info(f"Converting {converter}")
    start = datetime.now()
    try:
        update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.CONVERTING})
        file_dir = converter(local_filename)
        update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.CONVERTED})
        logger.info(f"Finished converting {converter} in {datetime.now()- start}")
    except Exception:
        file_dir = None
        update_db(dataset_id, processing_status={processing_status_type: ConversionStatus.FAILED})
        logger.exception(error_message)
    return file_dir


def get_bucket_prefix(dataset_id):
    remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
    if remote_dev_prefix:
        return join(remote_dev_prefix, dataset_id).strip("/")
    else:
        return dataset_id


def process_cxg(local_filename, dataset_id, cellxgene_bucket):
    bucket_prefix = get_bucket_prefix(dataset_id)
    cxg_dir = convert_file_ignore_exceptions(make_cxg, local_filename, "Issue creating cxg.", dataset_id, "cxg_status")
    if cxg_dir:
        update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.UPLOADING})
        s3_uri = copy_cxg_files_to_cxg_bucket(cxg_dir, bucket_prefix, cellxgene_bucket)
        metadata = {
            "explorer_url": join(DEPLOYMENT_STAGE_TO_URL[os.environ["DEPLOYMENT_STAGE"]], dataset_id + ".cxg", "")
        }
        with db_session_manager() as session:
            logger.info(f"Updating database with cxg artifact for dataset {dataset_id}. s3_uri is {s3_uri}")
            DatasetAsset.create(
                session,
                dataset_id=dataset_id,
                filename="explorer_cxg",
                filetype=DatasetArtifactFileType.CXG,
                type_enum=DatasetArtifactType.REMIX,
                user_submitted=True,
                s3_uri=s3_uri,
            )
        update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.UPLOADED})

    else:
        metadata = None
    update_db(dataset_id, metadata)


def validate_h5ad_file_and_add_labels(dataset_id, local_filename):
    from cellxgene_schema import validate

    update_db(dataset_id, processing_status=dict(validation_status=ValidationStatus.VALIDATING))
    output_filename = "local.h5ad"
    is_valid, errors = validate.validate(local_filename, output_filename)

    if not is_valid:
        logger.error(f"Validation failed with {len(errors)} errors!")
        status = dict(
            validation_status=ValidationStatus.INVALID,
            validation_message=errors,
            processing_status=ProcessingStatus.FAILURE,
        )
        update_db(dataset_id, processing_status=status)
        raise ValidationFailed
    else:
        logger.info("Validation complete")
        status = dict(
            h5ad_status=ConversionStatus.CONVERTED,
            validation_status=ValidationStatus.VALID,
        )
        update_db(dataset_id, processing_status=status)
        return output_filename


def clean_up_local_file(local_filename):
    try:
        os.remove(local_filename)
    except Exception:
        pass


def log_batch_environment():
    batch_environment_variables = [
        "AWS_BATCH_CE_NAME",
        "AWS_BATCH_JOB_ATTEMPT",
        "AWS_BATCH_JOB_ID",
        "DROPBOX_URL",
        "ARTIFACT_BUCKET",
        "CELLXGENE_BUCKET",
        "DATASET_ID",
        "DEPLOYMENT_STAGE",
    ]
    env_vars = dict()
    for var in batch_environment_variables:
        env_vars[var] = os.getenv(var)
    logger.info(f"Batch Job Info: {env_vars}")


def process(dataset_id, dropbox_url, cellxgene_bucket, artifact_bucket):
    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.PENDING))
    local_filename = download_from_dropbox_url(
        dataset_id,
        dropbox_url,
        "raw.h5ad",
    )

    # No file cleanup needed due to docker run-time environment.
    # To implement proper cleanup, tests/unit/backend/corpora/dataset_processing/test_process.py
    # will have to be modified since it relies on a shared local file

    file_with_labels = validate_h5ad_file_and_add_labels(dataset_id, local_filename)
    # Process metadata
    metadata = extract_metadata(file_with_labels)
    update_db(dataset_id, metadata)

    # create artifacts
    process_cxg(file_with_labels, dataset_id, cellxgene_bucket)
    create_artifacts(file_with_labels, dataset_id, artifact_bucket)
    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.SUCCESS))


def main():
    return_value = 0
    check_env()
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    try:
        process(dataset_id, os.environ["DROPBOX_URL"], os.environ["CELLXGENE_BUCKET"], os.environ["ARTIFACT_BUCKET"])
    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed):
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception:
        message = "An unexpected error occurred while processing the data set."
        logger.exception(message)
        update_db(
            dataset_id, processing_status=dict(processing_status=ProcessingStatus.FAILURE, upload_message=message)
        )
        return_value = 1

    if return_value > 0:
        notify_slack_failure(dataset_id)
    return return_value


def notify_slack_failure(dataset_id):
    data = format_slack_message(dataset_id)
    logger.info(data)
    if os.getenv("DEPLOYMENT_STAGE") == "prod":
        slack_webhook = CorporaConfig().slack_webhook
        requests.post(slack_webhook, headers={"Content-type": "application/json"}, data=data)


if __name__ == "__main__":
    rv = main()
    sys.exit(rv)
