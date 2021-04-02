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
    conversion_loom_status = ConversionStatus
    conversion_rds_status = ConversionStatus
    conversion_cxg_status = ConversionStatus
    conversion_anndata_status = ConversionStatus
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
    conversion_loom_status = ConversionStatus.CONVERTING
    conversion_rds_status = ConversionStatus.CONVERTING
    conversion_cxg_status = ConversionStatus.CONVERTING
    conversion_anndata_status = ConversionStatus.CONVERTING
}

If validation fails the processing_status change to:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.FAILED
}

## Conversion
After each conversion the processing_status change from CONVERTING to CONVERTED. Cellxgene data is converted first.
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    conversion_loom_status = ConversionStatus.CONVERTING
    conversion_rds_status = ConversionStatus.CONVERTING
    conversion_cxg_status = ConversionStatus.CONVERTED
    conversion_anndata_status = ConversionStatus.CONVERTING
}

If a conversion fails the processing_status will indicated it as follow:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    conversion_loom_status = ConversionStatus.FAILED
    conversion_rds_status = ConversionStatus.CONVERTING
    conversion_cxg_status = ConversionStatus.CONVERTED
    conversion_anndata_status = ConversionStatus.CONVERTING
}

Once all conversion are compelete, the conversion status for each file will be either CONVERTED or FAILED:
{
    processing_status = ProcessingStatus.SUCCESS
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    conversion_loom_status = ConversionStatus.FAILED
    conversion_rds_status = ConversionStatus.CONVERTED
    conversion_cxg_status = ConversionStatus.CONVERTED
    conversion_anndata_status = ConversionStatus.FAILED
}
"""

import logging
import os
import subprocess
import typing
from os.path import basename, join

import boto3
import numpy
import scanpy
import sys

from backend.corpora.common.utils.dropbox import get_download_url_from_shared_link, get_file_info
from backend.corpora.dataset_processing.exceptions import ProcessingFailed, ValidationFailed, ProcessingCancelled
from backend.corpora.common.corpora_orm import (
    DatasetArtifactFileType,
    DatasetArtifactType,
    ConversionStatus,
    ValidationStatus,
    ProcessingStatus,
)
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager, processing_status_updater
from backend.corpora.dataset_processing.download import download

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# This is unfortunate, but this information doesn't appear to live anywhere
# accessible to the uploader
DEPLOYMENT_STAGE_TO_URL = {
    "dev": "https://cellxgene.dev.single-cell.czi.technology/e",
    "staging": "https://cellxgene.staging.single-cell.czi.technology/e",
    "prod": "https://cellxgene.cziscience.com/e",
    "rdev": os.environ.get("FRONTEND_URL"),
}

s3_client = boto3.client(
    "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None, config=boto3.session.Config(signature_version="s3v4")
)


def check_env():
    """Verify that the required environment variables are set."""

    missing = []
    for env_var in ["DROPBOX_URL", "ARTIFACT_BUCKET", "CELLXGENE_BUCKET", "DATASET_ID", "DEPLOYMENT_STAGE"]:
        if env_var not in os.environ:
            missing.append(env_var)
    if missing:
        raise EnvironmentError(f"Missing environment variables: {missing}")


def create_artifact(
    file_name: str, artifact_type: DatasetArtifactFileType, bucket_prefix: str, dataset_id: str, artifact_bucket: str
) -> DatasetAsset:
    file_base = basename(file_name)
    logger.info(f"Uploading to [{artifact_type}] to S3 bucket: [{artifact_bucket}].")
    s3_client.upload_file(
        file_name,
        artifact_bucket,
        join(bucket_prefix, file_base),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )
    logger.info(f"Updating database with  {artifact_type}.")
    with db_session_manager() as session:
        DatasetAsset.create(
            session,
            dataset_id=dataset_id,
            filename=file_base,
            filetype=artifact_type,
            type_enum=DatasetArtifactType.REMIX,
            user_submitted=True,
            s3_uri=join("s3://", artifact_bucket, bucket_prefix, file_base),
        )


def create_artifacts(local_filename, dataset_id, artifact_bucket):
    bucket_prefix = get_bucket_prefix(dataset_id)
    logger.info(f"Creating Artifacts.")
    # upload AnnData
    create_artifact(local_filename, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket)
    update_db(
        dataset_id,
        processing_status=dict(conversion_anndata_status=ConversionStatus.CONVERTED),
    )

    # Process loom
    loom_filename, status = convert_file_ignore_exceptions(make_loom, local_filename, "Issue creating loom.")
    if loom_filename:
        create_artifact(loom_filename, DatasetArtifactFileType.LOOM, bucket_prefix, dataset_id, artifact_bucket)
    update_db(
        dataset_id,
        processing_status=dict(conversion_loom_status=status),
    )

    # Process seurat
    seurat_filename, status = convert_file_ignore_exceptions(make_seurat, local_filename, "Issue creating seurat.")
    if seurat_filename:
        create_artifact(seurat_filename, DatasetArtifactFileType.RDS, bucket_prefix, dataset_id, artifact_bucket)
    update_db(dataset_id, processing_status=dict(conversion_rds_status=status))


def cancel_dataset(dataset_id):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        dataset.dataset_and_asset_deletion()
        logger.info("Upload Canceled.")


def update_db(dataset_id, metadata=None, processing_status=None):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        if dataset.tombstone:
            raise ProcessingCancelled

        if metadata:
            # TODO: Delete this line once mean_genes_per_cell is in the db
            metadata.pop("mean_genes_per_cell", None)
            logger.debug(f"updating metadata.")
            dataset.update(**metadata)

        if processing_status:
            logger.debug(f"updating processing_status.{processing_status}")
            processing_status_updater(session, dataset.processing_status.id, processing_status)


def download_from_dropbox_url(dataset_uuid: str, dropbox_url: str, local_path: str) -> str:
    """Given a dropbox url, download it to local_path.
    Handles fixing the url so it downloads directly.
    """

    fixed_dropbox_url = get_download_url_from_shared_link(dropbox_url)
    if not fixed_dropbox_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    file_info = get_file_info(fixed_dropbox_url)
    status = download(dataset_uuid, fixed_dropbox_url, local_path, file_info["size"])
    logger.info(status)
    return local_path


def extract_metadata(filename):
    """Pull metadata out of the AnnData file to insert into the dataset table."""

    adata = scanpy.read_h5ad(filename, backed="r")

    try:
        raw_layer_name = [k for k, v in adata.uns["layer_descriptions"].items() if v == "raw"][0]
    except (KeyError, IndexError):
        raise RuntimeError("Raw layer not found in layer descriptions!")

    if raw_layer_name == "X":
        raw_layer = adata.X
    elif raw_layer_name == "raw.X":
        raw_layer = adata.raw.X
    else:
        raw_layer = adata.layers[raw_layer_name]

    # Calling np.count_nonzero on and h5py.Dataset appears to read the entire thing
    # into memory, so we need to chunk it to be safe.
    stride = 50000
    numerator, denominator = 0, 0
    for bounds in zip(range(0, raw_layer.shape[0], stride), range(stride, raw_layer.shape[0] + stride, stride)):
        chunk = raw_layer[bounds[0] : bounds[1], :]
        numerator += chunk.nnz if hasattr(chunk, "nnz") else numpy.count_nonzero(chunk)
        denominator += chunk.shape[0]

    def _get_term_pairs(base_term):
        base_term_id = base_term + "_ontology_term_id"
        return [
            {"label": k[0], "ontology_term_id": k[1]}
            for k in adata.obs.groupby([base_term, base_term_id]).groups.keys()
        ]

    return {
        "name": adata.uns["title"],
        "organism": {"label": adata.uns["organism"], "ontology_term_id": adata.uns["organism_ontology_term_id"]},
        "tissue": _get_term_pairs("tissue"),
        "assay": _get_term_pairs("assay"),
        "disease": _get_term_pairs("disease"),
        "sex": list(adata.obs.sex.unique()),
        "ethnicity": _get_term_pairs("ethnicity"),
        "development_stage": _get_term_pairs("development_stage"),
        "cell_count": adata.shape[0],
        "mean_genes_per_cell": numerator / denominator,
    }


def make_loom(local_filename):
    """Create a loom file from the AnnData file."""

    adata = scanpy.read_h5ad(local_filename)
    column_name_map = {}
    for column in adata.obs.columns:
        if "/" in column:
            column_name_map[column] = column.replace("/", "-")
    if column_name_map:
        adata.obs = adata.obs.rename(columns=column_name_map)

    loom_filename = local_filename.replace(".h5ad", ".loom")
    adata.write_loom(loom_filename, True)
    return loom_filename


def make_seurat(local_filename):
    """Create a Seurat rds file from the AnnData file."""

    try:
        seurat_proc = subprocess.run(
            ["Rscript", os.path.join(os.path.abspath(os.path.dirname(__file__)), "make_seurat.R"), local_filename],
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        logger.exception()
        raise RuntimeError(f"Seurat conversion failed: {seurat_proc.stdout} {seurat_proc.stderr}")

    return local_filename.replace(".h5ad", ".rds")


def make_cxg(local_filename):
    cxg_dir = local_filename.replace(".h5ad", ".cxg")
    cxg_proc = subprocess.run(
        ["cellxgene", "convert", "-o", cxg_dir, "-s", "10.0", local_filename], capture_output=True
    )
    if cxg_proc.returncode != 0:
        raise RuntimeError(f"CXG conversion failed: {cxg_proc.stderr}")
    return cxg_dir


def copy_cxg_files_to_cxg_bucket(cxg_dir, bucket_prefix, cellxgene_bucket):
    command = ["aws"]
    if os.getenv("BOTO_ENDPOINT_URL"):
        command.append(f"--endpoint-url={os.getenv('BOTO_ENDPOINT_URL')}")
    command.extend(
        [
            "s3",
            "cp",
            cxg_dir,
            f"s3://{cellxgene_bucket}/{bucket_prefix}.cxg/",
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
    converter: typing.Callable, local_filename: str, error_message: str
) -> typing.Tuple[str, ConversionStatus]:
    try:
        file_dir = converter(local_filename)
        status = ConversionStatus.CONVERTED
    except Exception:
        file_dir = None
        status = ConversionStatus.FAILED
        logger.exception(error_message)
    return file_dir, status


def get_bucket_prefix(dataset_id):
    remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
    if remote_dev_prefix:
        return join(remote_dev_prefix, dataset_id).strip("/")
    else:
        return dataset_id


def process_cxg(local_filename, dataset_id, cellxgene_bucket):
    bucket_prefix = get_bucket_prefix(dataset_id)
    cxg_dir, status = convert_file_ignore_exceptions(make_cxg, local_filename, "Issue creating cxg.")
    if cxg_dir:
        copy_cxg_files_to_cxg_bucket(cxg_dir, bucket_prefix, cellxgene_bucket)
        metadata = {
            "deployment_directories": [
                {"url": join(DEPLOYMENT_STAGE_TO_URL[os.environ["DEPLOYMENT_STAGE"]], dataset_id + ".cxg", "")}
            ]
        }
    else:
        metadata = None
    update_db(dataset_id, metadata, processing_status=dict(conversion_cxg_status=status))


def validate_h5ad_file(dataset_id, local_filename):
    update_db(dataset_id, processing_status=dict(validation_status=ValidationStatus.VALIDATING))
    val_proc = subprocess.run(["cellxgene-schema", "validate", local_filename], capture_output=True)
    if val_proc.returncode != 0:
        logger.error("Validation failed!")
        logger.error(f"stdout: {val_proc.stdout}")
        logger.error(f"stderr: {val_proc.stderr}")
        status = dict(
            validation_status=ValidationStatus.INVALID,
            validation_message=val_proc.stdout,
            processing_status=ProcessingStatus.FAILURE,
        )
        update_db(dataset_id, processing_status=status)
        raise ValidationFailed
    else:
        logger.info("Validation complete", flush=True)
        status = dict(
            conversion_cxg_status=ConversionStatus.CONVERTING,
            conversion_loom_status=ConversionStatus.CONVERTING,
            conversion_rds_status=ConversionStatus.CONVERTING,
            conversion_anndata_status=ConversionStatus.CONVERTING,
            validation_status=ValidationStatus.VALID,
        )
        update_db(dataset_id, processing_status=status)


def log_batch_environment():
    batch_environment_variables = ["AWS_BATCH_CE_NAME", "AWS_BATCH_JOB_ATTEMPT", "AWS_BATCH_JOB_ID"]
    vars = dict()
    for var in batch_environment_variables:
        vars[var] = os.getenv(var)
    logger.info(f"Batch Job Info: {vars}")


def main():
    check_env()
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    try:
        update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.PENDING))
        local_filename = download_from_dropbox_url(
            dataset_id,
            os.environ["DROPBOX_URL"],
            "local.h5ad",
        )
        logger.info("Download complete", flush=True)

        validate_h5ad_file(dataset_id, local_filename)
        process_cxg(local_filename, dataset_id, os.environ["CELLXGENE_BUCKET"])

        # Process metadata
        metadata = extract_metadata(local_filename)
        logger.info(metadata, flush=True)
        update_db(dataset_id, metadata)

        # create artifacts
        create_artifacts(
            local_filename,
            dataset_id,
            os.environ["ARTIFACT_BUCKET"],
        )
        update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.SUCCESS))
    except ProcessingCancelled:
        cancel_dataset(dataset_id)
        sys.exit(0)
    except ProcessingFailed:
        logging.exception()
        update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.FAILURE))
        sys.exit(1)
    except ValidationFailed:
        logger.exception()
        sys.exit(1)
    except Exception:
        message = "An unexpect error occured while processing the data set."
        logger.exception(message)
        update_db(
            dataset_id, processing_status=dict(processing_status=ProcessingStatus.FAILURE, upload_message=message)
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
