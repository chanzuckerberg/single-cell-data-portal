#!/usr/bin/env python3
import os
import subprocess
from os.path import join

from backend.common.corpora_orm import DatasetArtifactFileType, ConversionStatus
from backend.common.entities import DatasetAsset
from backend.common.utils.db_session import db_session_manager
from backend.dataset_processing.h5ad_data_file import H5ADDataFile
from backend.dataset_processing.logger import logger
from backend.dataset_processing.common import (
    download_from_s3,
    get_bucket_prefix,
    update_db,
    convert_file,
    DEPLOYMENT_STAGE_TO_URL,
)


def process(dataset_id: str, artifact_bucket: str, cellxgene_bucket: str):
    """
    1. Download the labeled dataset from the artifact bucket
    2. Convert the labeled dataset to CXG
    3. Upload the CXG to the cellxgene bucket
    :param dataset_id:
    :param artifact_bucket:
    :param cellxgene_bucket:
    :return:
    """

    labeled_h5ad_filename = "local.h5ad"

    # Download the labeled dataset from the artifact bucket
    bucket_prefix = get_bucket_prefix(dataset_id)
    object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
    download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

    # Convert the labeled dataset to CXG and upload it to the cellxgene bucket
    process_cxg(labeled_h5ad_filename, dataset_id, cellxgene_bucket)


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


def process_cxg(local_filename, dataset_id, cellxgene_bucket):
    cxg_dir = convert_file(make_cxg, local_filename, "Issue creating cxg.", dataset_id, "cxg_status")
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
