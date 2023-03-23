import logging
import os
import subprocess
from typing import Dict, Union

import anndata
import numpy as np
import requests
from anndata._core.views import ArrayView
from scipy import sparse

from backend.wmg.data.constants import INCLUDED_ASSAYS

logger = logging.getLogger(__name__)


def get_X_raw(anndata_object: anndata.AnnData) -> Union[np.ndarray, sparse.spmatrix, ArrayView]:
    """
    Current rules for our curated H5ADs:
    * if there is a .raw, it is the raw counts, and .X is transformed/normalized (by author) or is == to .raw.X
    * if there is no .raw, ad.X contains the raw counts
    """
    raw_expression_matrix = getattr(anndata_object.raw, "X", None)
    return raw_expression_matrix if raw_expression_matrix is not None else anndata_object.X


def get_dataset_s3_uris() -> Dict[str, str]:
    """
    Retrieve list of s3 uris for datasets included in the wmg cube
    """
    # hardcode to dev backend if deployment is rdev
    API_URL = (
        "https://api.cellxgene.dev.single-cell.czi.technology"
        if os.environ.get("DEPLOYMENT_STAGE") in ["test", "rdev"]
        else os.getenv("API_URL")
    )

    dataset_metadata_url = f"{API_URL}/dp/v1/datasets/index"
    datasets = requests.get(dataset_metadata_url).json()

    s3_uris = dict()
    for dataset in datasets:
        if (
            dataset["is_primary_data"] == "PRIMARY"
            and any(assay["ontology_term_id"] in INCLUDED_ASSAYS for assay in dataset["assay"])
            and len(dataset["organism"]) < 2
        ):
            dataset_id = dataset["explorer_url"].split("/")[-2].split(".cxg")[0]
            s3_uri = next(
                a["s3_uri"]
                for a in dataset["dataset_assets"]
                if a["filetype"] == "H5AD" and a["filename"] == "local.h5ad"
            )
            s3_uris[dataset_id] = s3_uri
    return s3_uris


def copy_datasets_to_instance(s3_uris: Dict, dataset_directory: str):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        copy_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
        subprocess.run(copy_command)


def extract_h5ad(h5ad_path: str) -> anndata.AnnData:
    logger.info(f"Extracting {h5ad_path}...")
    return anndata.read_h5ad(h5ad_path)
