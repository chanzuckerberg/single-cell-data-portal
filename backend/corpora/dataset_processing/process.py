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

# Standalone processing steps

## Seurat
This is used to recompute the Seurat artifact in place, starting from the original h5ad.
This is a state machine with a single state that mimics the Conversion step
of the main step function.

## CXG_Remaster
This is used to migrate the cxg to a different, more performant format. This is a state machine with a single
state that mimics the Conversion step of the main step function.

"""

import logging
import os
import subprocess
import sys
import typing
from datetime import datetime
from os.path import join

import numpy
import scanpy

from backend.corpora.common.corpora_orm import (
    ConversionStatus,
    DatasetArtifactFileType,
    ProcessingStatus,
    ValidationStatus,
    UploadStatus,
)
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.dl_sources.url import from_url
from backend.corpora.common.utils.s3_buckets import buckets
from backend.corpora.dataset_processing.download import download
from backend.corpora.dataset_processing.exceptions import (
    ProcessingCancelled,
    ProcessingFailed,
    ValidationFailed,
)
from backend.corpora.dataset_processing.h5ad_data_file import H5ADDataFile

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# This is unfortunate, but this information doesn't appear to live anywhere
# accessible to the uploader
DEPLOYMENT_STAGE_TO_URL = {
    "staging": "https://cellxgene.staging.single-cell.czi.technology/e",
    "prod": "https://cellxgene.cziscience.com/e",
    "rdev": f"https:/{os.environ.get('REMOTE_DEV_PREFIX')}-explorer.rdev.single-cell.czi.technology/e",
    "dev": "https://cellxgene.dev.single-cell.czi.technology/e",
    "test": "http://frontend.corporanet.local:3000",
}

LABELED_H5AD_FILENAME = "local.h5ad"


def check_env():
    """Verify that the required environment variables are set."""

    missing = []
    for env_var in [
        "DROPBOX_URL",
        "ARTIFACT_BUCKET",
        "CELLXGENE_BUCKET",
        "DATASET_ID",
        "DEPLOYMENT_STAGE",
    ]:
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
    update_db(
        dataset_id,
        processing_status={processing_status_type: ConversionStatus.UPLOADING},
    )
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
                user_submitted=True,
                s3_uri=s3_uri,
            )
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.UPLOADED},
        )

    except Exception as e:
        logger.error(e)
        e.args = {processing_status_type: ConversionStatus.FAILED}
        raise e


def replace_artifact(
    file_name: str,
    bucket_prefix: str,
    artifact_bucket: str,
):
    logger.info(f"Uploading [{bucket_prefix}/{file_name}] to S3 bucket: [{artifact_bucket}].")
    DatasetAsset.upload(file_name, bucket_prefix, artifact_bucket)


def create_artifacts(
    local_filename: str,
    dataset_id: str,
    artifact_bucket: str,
    can_convert_to_seurat: bool = False,
):
    bucket_prefix = get_bucket_prefix(dataset_id)
    logger.info(f"Creating artifacts for dataset {dataset_id}...")
    # upload AnnData
    create_artifact(
        local_filename,
        DatasetArtifactFileType.H5AD,
        bucket_prefix,
        dataset_id,
        artifact_bucket,
        "h5ad_status",
    )

    if can_convert_to_seurat:
        # Convert to Seurat and upload
        seurat_filename = convert_file_ignore_exceptions(
            make_seurat,
            local_filename,
            "Issue creating seurat.",
            dataset_id,
            "rds_status",
        )
        if seurat_filename:
            create_artifact(
                seurat_filename,
                DatasetArtifactFileType.RDS,
                bucket_prefix,
                dataset_id,
                artifact_bucket,
                "rds_status",
            )
    else:
        update_db(dataset_id, processing_status=dict(rds_status=ConversionStatus.SKIPPED))
        logger.info(f"Skipped Seurat conversion for dataset {dataset_id}")

    logger.info(f"Finished creating artifacts for dataset {dataset_id}")


def cancel_dataset(dataset_id):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        dataset.asset_deletion()
        dataset.delete()
        logger.info("Upload Canceled.")


def update_db(dataset_id, metadata: dict = None, processing_status: dict = None):
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


def download_from_source_uri(dataset_id: str, source_uri: str, local_path: str) -> str:
    """Given a source URI, download it to local_path.
    Handles fixing the url so it downloads directly.
    """

    file_url = from_url(source_uri)
    if not file_url:
        raise ValueError(f"Malformed source URI: {source_uri}")

    # This is a bit ugly and should be done polymorphically instead, but Dropbox support will be dropped soon
    if file_url.scheme == "https":
        file_info = file_url.file_info()
        status = download(dataset_id, file_url.url, local_path, file_info["size"])
        logger.info(status)
    elif file_url.scheme == "s3":
        bucket_name = file_url.netloc
        key = remove_prefix(file_url.path, "/")
        wrapped_download_from_s3(
            dataset_id=dataset_id,
            bucket_name=bucket_name,
            object_key=key,
            local_filename=local_path,
        )
    else:
        raise ValueError(f"Download for URI scheme '{file_url.scheme}' not implemented")

    logger.info("Download complete")
    return local_path


# TODO: after upgrading to Python 3.9, replace this with removeprefix()
def remove_prefix(string: str, prefix: str) -> str:
    if string.startswith(prefix):
        return string[len(prefix) :]
    else:
        return string[:]


def download_from_s3(bucket_name: str, object_key: str, local_filename: str):
    logger.info(f"Downloading file {local_filename} from bucket {bucket_name} with object key {object_key}")
    buckets.portal_client.download_file(bucket_name, object_key, local_filename)


def wrapped_download_from_s3(dataset_id: str, bucket_name: str, object_key: str, local_filename: str):
    """
    Wraps download_from_s3() to update the dataset's upload status
    :param dataset_id:
    :param bucket_name:
    :param object_key:
    :param local_filename:
    :return:
    """
    with db_session_manager() as session:
        processing_status = Dataset.get(session, dataset_id).processing_status
        processing_status.upload_status = UploadStatus.UPLOADING
        download_from_s3(
            bucket_name=bucket_name,
            object_key=object_key,
            local_filename=local_filename,
        )
        processing_status.upload_status = UploadStatus.UPLOADED


def extract_metadata(filename) -> dict:
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

    def _get_batch_condition():
        if "batch_condition" in adata.uns:
            return adata.uns["batch_condition"]
        else:
            return None

    metadata = {
        "name": adata.uns["title"],
        "organism": _get_term_pairs("organism"),
        "tissue": _get_term_pairs("tissue"),
        "assay": _get_term_pairs("assay"),
        "disease": _get_term_pairs("disease"),
        "sex": _get_term_pairs("sex"),
        "ethnicity": _get_term_pairs("self_reported_ethnicity"),
        "development_stage": _get_term_pairs("development_stage"),
        "cell_count": adata.shape[0],
        "mean_genes_per_cell": numerator / denominator,
        "is_primary_data": _get_is_primary_data(),
        "cell_type": _get_term_pairs("cell_type"),
        "x_approximate_distribution": _get_x_approximate_distribution(),
        "schema_version": adata.uns["schema_version"],
        "batch_condition": _get_batch_condition(),
        "donor_id": adata.obs["donor_id"].unique(),
        "suspension_type": adata.obs["suspension_type"].unique(),
    }
    logger.info(f"Extract metadata: {metadata}")
    return metadata


def make_seurat(local_filename):
    """Create a Seurat rds file from the AnnData file."""

    try:
        subprocess.run(
            [
                "Rscript",
                os.path.join(os.path.abspath(os.path.dirname(__file__)), "make_seurat.R"),
                local_filename,
            ],
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
        h5ad_data_file.to_cxg(cxg_output_container, sparse_threshold=25.0)
    except Exception as ex:
        msg = "CXG conversion failed."
        logger.exception(msg)
        raise RuntimeError(msg) from ex

    return cxg_output_container


def copy_cxg_files_to_cxg_bucket(cxg_dir, s3_uri):
    """
    Copy cxg files to the cellxgene bucket (under the given object key) for access by the explorer
    """
    command = ["aws"]
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


def convert_file_ignore_exceptions(
    converter: typing.Callable,
    local_filename: str,
    error_message: str,
    dataset_id: str,
    processing_status_type: str,
) -> str:
    logger.info(f"Converting {local_filename}")
    start = datetime.now()
    try:
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.CONVERTING},
        )
        file_dir = converter(local_filename)
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.CONVERTED},
        )
        logger.info(f"Finished converting {converter} in {datetime.now()- start}")
    except Exception:
        file_dir = None
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.FAILED},
        )
        logger.exception(error_message)
    return file_dir


def get_bucket_prefix(identifier):
    remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
    if remote_dev_prefix:
        return join(remote_dev_prefix, identifier).strip("/")
    else:
        return identifier


def process_cxg(local_filename, dataset_id, cellxgene_bucket):
    cxg_dir = convert_file_ignore_exceptions(make_cxg, local_filename, "Issue creating cxg.", dataset_id, "cxg_status")
    if cxg_dir:
        with db_session_manager() as session:
            asset = DatasetAsset.create(
                session,
                dataset_id=dataset_id,
                filename="explorer_cxg",
                filetype=DatasetArtifactFileType.CXG,
                user_submitted=True,
                s3_uri="",
            )
            asset_id = asset.id
            bucket_prefix = get_bucket_prefix(asset_id)
            s3_uri = f"s3://{cellxgene_bucket}/{bucket_prefix}.cxg/"
        update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.UPLOADING})
        copy_cxg_files_to_cxg_bucket(cxg_dir, s3_uri)
        metadata = {
            "explorer_url": join(
                DEPLOYMENT_STAGE_TO_URL[os.environ["DEPLOYMENT_STAGE"]],
                dataset_id + ".cxg",
                "",
            )
        }
        with db_session_manager() as session:
            logger.info(f"Updating database with cxg artifact for dataset {dataset_id}. s3_uri is {s3_uri}")
            asset = DatasetAsset.get(session, asset_id)
            asset.update(s3_uri=s3_uri)
        update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.UPLOADED})

    else:
        metadata = None
    update_db(dataset_id, metadata)


def validate_h5ad_file_and_add_labels(dataset_id: str, local_filename: str) -> typing.Tuple[str, bool]:
    """
    Validates and labels the specified dataset file and updates the processing status in the database
    :param dataset_id: ID of the dataset to update
    :param local_filename: file name of the dataset to validate and label
    :return: file name of labeled dataset, boolean indicating if seurat conversion is possible
    """
    from cellxgene_schema import validate

    update_db(
        dataset_id,
        processing_status=dict(validation_status=ValidationStatus.VALIDATING),
    )
    output_filename = LABELED_H5AD_FILENAME
    is_valid, errors, can_convert_to_seurat = validate.validate(local_filename, output_filename)

    if not is_valid:
        logger.error(f"Validation failed with {len(errors)} errors!")
        status = dict(
            validation_status=ValidationStatus.INVALID,
            validation_message=errors,
        )
        raise ValidationFailed(status)
    else:
        logger.info("Validation complete")
        status = dict(
            h5ad_status=ConversionStatus.CONVERTED,
            validation_status=ValidationStatus.VALID,
        )
        update_db(dataset_id, processing_status=status)
        return output_filename, can_convert_to_seurat


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
    logger.info(f"Batch Job Info: {env_vars}")


def process(dataset_id, dropbox_url, cellxgene_bucket, artifact_bucket):
    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.PENDING))
    local_filename = download_from_source_uri(
        dataset_id=dataset_id,
        source_uri=dropbox_url,
        local_path="raw.h5ad",
    )

    # No file cleanup needed due to docker run-time environment.
    # To implement proper cleanup, tests/unit/backend/corpora/dataset_processing/test_process.py
    # will have to be modified since it relies on a shared local file

    file_with_labels, can_convert_to_seurat = validate_h5ad_file_and_add_labels(dataset_id, local_filename)

    # Process metadata
    metadata = extract_metadata(file_with_labels)
    update_db(dataset_id, metadata)

    # create artifacts
    process_cxg(file_with_labels, dataset_id, cellxgene_bucket)
    create_artifacts(file_with_labels, dataset_id, artifact_bucket, can_convert_to_seurat)
    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.SUCCESS))


def main():
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    step_name = os.environ["STEP_NAME"]
    is_last_attempt = os.environ["AWS_BATCH_JOB_ATTEMPT"] == os.getenv("MAX_ATTEMPTS", "1")
    return_value = 0
    logger.info(f"Processing dataset {dataset_id}")
    try:
        if step_name == "download-validate":
            from backend.corpora.dataset_processing.process_download_validate import (
                process,
            )

            process(dataset_id, os.environ["DROPBOX_URL"], os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg":
            from backend.corpora.dataset_processing.process_cxg import process

            process(
                dataset_id,
                os.environ["ARTIFACT_BUCKET"],
                os.environ["CELLXGENE_BUCKET"],
            )
        elif step_name == "seurat":
            from backend.corpora.dataset_processing.process_seurat import process

            process(dataset_id, os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg_remaster":
            from backend.corpora.dataset_processing.remaster_cxg import process

            process(dataset_id, os.environ["CELLXGENE_BUCKET"], dry_run=False)
        else:
            logger.error(f"Step function configuration error: Unexpected STEP_NAME '{step_name}'")

    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed) as e:
        (status,) = e.args
        if is_last_attempt:
            update_db(dataset_id, processing_status=status)
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception as e:
        (status,) = e.args
        if is_last_attempt and isinstance(status, dict):
            update_db(dataset_id, processing_status=status)
        logger.exception(f"An unexpected error occurred while processing the data set: {e}")
        return_value = 1

    return return_value


if __name__ == "__main__":
    rv = main()
    sys.exit(rv)
