#!/usr/bin/env python3
import logging

from backend.corpora.common.utils.s3_buckets import s3_client
from backend.corpora.dataset_processing.process import (
    process_cxg,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def download_from_s3(bucket_name: str, object_key: str, local_filename: str):
    logger.info(f"Downloading file {local_filename} from bucket {bucket_name} with object key {object_key}")
    s3_client.download_file(bucket_name, object_key, local_filename)


def process(dataset_id: str, artifact_bucket: str, cellxgene_bucket: str):
    """
    1. Download the labeled dataset from the artifact bucket
    2. Convert the labeled dataset to CXG
    3. Upload the CXG to the cellxgene bucket
    :param dataset_id:
    :param artifact_bucket:
    :param labeled_h5ad_filename:
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
