import logging
import subprocess
from typing import Dict, Union

import anndata
import numpy as np
from anndata._core.views import ArrayView
from scipy import sparse

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetArtifactType
from backend.layers.persistence.persistence import DatabaseProvider
from backend.wmg.data.constants import INCLUDED_ASSAYS

logger = logging.getLogger(__name__)

_business_logic = None


def get_business_logic():
    """
    Returns an instance of the business logic handler. Use this to interrogate the database
    """
    global _business_logic
    if not _business_logic:
        _business_logic = BusinessLogic(DatabaseProvider(), None, None, None, None)
    return _business_logic


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
    s3_uris = dict()
    for dataset in get_business_logic().get_all_published_datasets():
        if (
            (dataset.metadata is not None)
            and (dataset.metadata.assay is not None)
            and (dataset.metadata.is_primary_data == "PRIMARY")
            and (dataset.metadata.organism is not None)
            and any(assay.ontology_term_id in INCLUDED_ASSAYS for assay in dataset.metadata.assay)
            and len(dataset.metadata.organism) < 2
        ):
            s3_uri = next(
                a.uri
                for a in dataset.artifacts
                if a.type == DatasetArtifactType.H5AD and a.get_file_name() == "local.h5ad"
            )
            s3_uris[dataset.dataset_id.id] = s3_uri
    return s3_uris


def copy_datasets_to_instance(s3_uris: Dict, dataset_directory: str):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        copy_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
        subprocess.run(copy_command)


def extract_h5ad(h5ad_path: str) -> anndata.AnnData:
    logger.info(f"Extracting {h5ad_path}...")
    return anndata.read_h5ad(h5ad_path)
