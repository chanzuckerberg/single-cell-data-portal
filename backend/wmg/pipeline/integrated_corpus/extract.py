import logging
import subprocess
from typing import Dict, Union

import anndata
import numpy as np
from anndata._core.views import ArrayView
from scipy import sparse

from backend.common.corpora_orm import DatasetArtifactFileType
from backend.common.entities import Dataset, Collection, DatasetAsset
from backend.common.utils.db_session import db_session_manager
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


def get_dataset_s3_uris() -> Dict:
    """
    Retrieve list of s3 uris for datasets included in the wmg cube
    """
    with db_session_manager() as session:

        dataset_ids = []
        published_dataset_non_null_assays = (
            session.query(Dataset.table.id, Dataset.table.assay, Dataset.table.organism)
            .join(Dataset.table.collection)
            .filter(
                Dataset.table.assay != "null",
                Dataset.table.published == "TRUE",
                Dataset.table.is_primary_data == "PRIMARY",
                Collection.table.visibility == "PUBLIC",
                Dataset.table.tombstone == "FALSE",
                Dataset.table.organism != "null",
            )
            .all()
        )

        for dataset_id, assays, organisms in published_dataset_non_null_assays:
            if any(assay["ontology_term_id"] in INCLUDED_ASSAYS for assay in assays):
                if len(organisms) < 2:
                    dataset_ids.append(dataset_id)

        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD, "local.h5ad")
    return s3_uris


def copy_datasets_to_instance(s3_uris: Dict, dataset_directory: str):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        copy_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
        subprocess.run(copy_command)


def extract_h5ad(h5ad_path: str) -> anndata.AnnData:
    logger.info(f"Extracting {h5ad_path}...")
    return anndata.read_h5ad(h5ad_path)
