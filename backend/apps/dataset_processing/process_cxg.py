#!/usr/bin/env python3
import logging

from backend.apps.dataset_processing.process import (
    download_from_s3,
    process_cxg,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
