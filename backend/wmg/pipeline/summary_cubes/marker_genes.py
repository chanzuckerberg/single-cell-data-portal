import logging

import pandas as pd
import tiledb
import numpy as np

from backend.wmg.data.schemas.marker_genes_cube_schema import marker_genes_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, MARKER_GENES_CUBE_NAME
from backend.wmg.data.utils import create_empty_cube, log_func_runtime
from backend.wmg.api.calculate_markers import get_markers

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def extract_cellcounts(corpus_path: str) -> pd.DataFrame:
    """
    get obs data from integrated corpus
    """
    return tiledb.open(f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}")


@log_func_runtime
def create_marker_genes_cube(corpus_path: str):
    """
    Create marker genes cube and write to disk
    """
    cell_counts = extract_cellcounts()
    tissues_celltypes = (
        cell_counts.query(attrs=["cell_type_ontology_term_id"], dims=["tissue_ontology_term_id"])
        .df[:]
        .groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"])
        .first()
    )
    uniq_tissues = tissues.index.levels[0]
    tissues = np.array(list(tissues_celltypes.index.get_level_values(0)))
    cell_types = np.array(list(tissues_celltypes.index.get_level_values(1)))
    df = []
    for tiss in uniq_tissues:
        tiss_celltypes = cell_types[tissues == tiss]
        context = {
            "tissue_ontology_term_ids": [tiss],
            "cell_type_ontology_term_ids": tiss_celltypes,
        }
        for ct in tiss_celltypes:
            logger.info("Calculating markers for tissue: %s, cell type: %s", tiss, ct)
            target = {
                "tissue_ontology_term_ids": [tiss],
                "cell_type_ontology_term_ids": [ct],
            }
            t_markers = get_markers(target, context, corpus_path=corpus_path, test="ttest", n_markers=None)
            b_markers = get_markers(target, context, corpus_path=corpus_path, test="binomtest", n_markers=None)
            for g in t_markers:
                b_markers[g].update(t_markers[g])
            df.append(pd.from_json(b_markers))
    df = pd.concat(df)

    uri = f"{corpus_path}/{MARKER_GENES_CUBE_NAME}"
    create_empty_cube(uri, marker_genes_schema)
    tiledb.from_pandas(uri, df, mode="append")

    logger.info(f"Marker genes cube created and stored at {uri}")
