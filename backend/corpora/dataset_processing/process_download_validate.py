#!/usr/bin/env python3
import typing

import numpy
import scanpy

from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.dl_sources.url import from_url
from backend.corpora.dataset_processing.download import download

from backend.corpora.dataset_processing.exceptions import ValidationFailed
from backend.corpora.dataset_processing.common import create_artifact, update_db, download_from_s3, get_bucket_prefix
from backend.corpora.dataset_processing.logger import logger

from backend.corpora.common.corpora_orm import (
    ConversionStatus,
    ProcessingStatus,
    DatasetArtifactFileType,
    ValidationStatus,
    UploadStatus,
)


def process(dataset_id: str, dropbox_url: str, artifact_bucket: str):
    """
    1. Download the original dataset from Dropbox
    2. Validate and label it
    3. Upload the labeled dataset to the artifact bucket
    :param dataset_id:
    :param dropbox_url:
    :param artifact_bucket:
    :return:
    """

    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.PENDING))

    # Download the original dataset from Dropbox
    local_filename = download_from_source_uri(
        dataset_id=dataset_id,
        source_uri=dropbox_url,
        local_path="raw.h5ad",
    )

    # Validate and label the dataset
    file_with_labels, can_convert_to_seurat = validate_h5ad_file_and_add_labels(dataset_id, local_filename)
    # Process metadata
    metadata = extract_metadata(file_with_labels)
    update_db(dataset_id, metadata)

    if not can_convert_to_seurat:
        update_db(dataset_id, processing_status=dict(rds_status=ConversionStatus.SKIPPED))
        logger.info(f"Skipping Seurat conversion for dataset {dataset_id}")

    # Upload the labeled dataset to the artifact bucket
    bucket_prefix = get_bucket_prefix(dataset_id)
    create_artifact(
        file_with_labels, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )
    create_artifact(
        local_filename, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )


LABELED_H5AD_FILENAME = "local.h5ad"


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
    try:
        is_valid, errors, can_convert_to_seurat = validate.validate(local_filename, output_filename)
    except Exception as e:
        logger.error(f"Validation failed with exception: {e}!")
        status = dict(
            validation_status=ValidationStatus.INVALID,
            validation_message=str(e),
        )
        raise ValidationFailed(status)

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
        "ethnicity": _get_term_pairs("ethnicity"),
        "development_stage": _get_term_pairs("development_stage"),
        "cell_count": adata.shape[0],
        "mean_genes_per_cell": numerator / denominator,
        "is_primary_data": _get_is_primary_data(),
        "cell_type": _get_term_pairs("cell_type"),
        "x_normalization": adata.uns["X_normalization"],
        "x_approximate_distribution": _get_x_approximate_distribution(),
        "schema_version": adata.uns["schema_version"],
        "batch_condition": _get_batch_condition(),
    }
    logger.info(f"Extract metadata: {metadata}")
    return metadata


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
