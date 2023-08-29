import logging
import pathlib
import subprocess
from typing import Dict, Union

import anndata
import numpy as np
from anndata._core.views import ArrayView
from scipy import sparse

from backend.wmg.data.constants import INCLUDED_ASSAYS
from backend.wmg.data.utils import get_datasets_from_curation_api

logger: logging.Logger = logging.getLogger(__name__)


def get_X_raw(anndata_object: anndata.AnnData) -> Union[np.ndarray, sparse.spmatrix, ArrayView]:
    """
    Current rules for our curated H5ADs:
    * if there is a .raw, it is the raw counts, and .X is transformed/normalized (by author) or is == to .raw.X
    * if there is no .raw, ad.X contains the raw counts
    """
    raw_expression_matrix = getattr(anndata_object.raw, "X", None)
    return raw_expression_matrix if raw_expression_matrix is not None else anndata_object.X


def get_dataset_asset_urls(datasets=None) -> Dict[str, str]:
    """
    Retrieve list of asset urls for datasets included in the wmg cube

    :param datasets: list of datasets to check, if None, will retrieve from API
        This parameter is used for tests.
    """
    if datasets is None:
        datasets = get_datasets_from_curation_api()

    asset_urls = dict()
    for dataset in datasets:
        if (
            dataset["organism"] is not None
            and dataset["assay"] is not None
            and len(dataset["is_primary_data"]) == 1
            and dataset["is_primary_data"][0]
            and any(assay["ontology_term_id"] in INCLUDED_ASSAYS for assay in dataset["assay"])
            and len(dataset["organism"]) < 2
        ):
            dataset_id = dataset["dataset_id"]
            asset_url = next(a["url"] for a in dataset["assets"] if a["filetype"] == "H5AD")
            asset_urls[dataset_id] = asset_url
    return asset_urls


def copy_datasets_to_instance(s3_uris: Dict, dataset_directory: str):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        pathlib.Path(f"./{dataset_directory}/{dataset}").mkdir(parents=True, exist_ok=True)
        copy_command = ["wget", s3_uris[dataset], "-O", f"./{dataset_directory}/{dataset}/local.h5ad"]
        subprocess.run(copy_command)


def extract_h5ad(h5ad_path: str) -> anndata.AnnData:
    logger.info(f"Extracting {h5ad_path}...")
    return anndata.read_h5ad(h5ad_path)
