import logging

import pandas as pd
import tiledb
import numpy as np
import gc
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
    cell_counts = extract_cellcounts(corpus_path)
    tissues_celltypes = (
        cell_counts.query(
            attrs=["cell_type_ontology_term_id"], dims=["organism_ontology_term_id", "tissue_ontology_term_id"]
        )
        .df[:]
        .groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id", "organism_ontology_term_id"])
        .first()
    )
    uniq_tissues = tissues_celltypes.index.levels[0]
    tissues = np.array(list(tissues_celltypes.index.get_level_values(0)))
    cell_types = np.array(list(tissues_celltypes.index.get_level_values(1)))
    organisms = np.array(list(tissues_celltypes.index.get_level_values(2)))

    uri = f"{corpus_path}/{MARKER_GENES_CUBE_NAME}"
    create_empty_cube(uri, marker_genes_schema)

    for tiss in uniq_tissues:
        tiss_celltypes = list(cell_types[tissues == tiss])
        tiss_organisms = list(organisms[tissues == tiss])
        for i, ct in enumerate(tiss_celltypes):
            organism = tiss_organisms[i]

            logger.info("Calculating markers for tissue: %s, cell type: %s, organism: %s", tiss, ct, organism)
            target = {
                "tissue_ontology_term_ids": [tiss],
                "cell_type_ontology_term_ids": [ct],
                "organism_ontology_term_id": organism,
            }
            context = {
                "tissue_ontology_term_ids": [tiss],
                "organism_ontology_term_id": organism,
            }
            t_markers = get_markers(target, context, corpus=corpus_path, test="ttest", n_markers=None)
            b_markers = get_markers(target, context, corpus=corpus_path, test="binomtest", n_markers=None)
            gc.collect()

            all_marker_genes = set(list(t_markers.keys())).union(list(b_markers.keys()))
            markers = []
            for g in all_marker_genes:
                b_stats = b_markers.get(g, {"p_value_binomtest": 1, "effect_size_binomtest": 0})
                t_stats = t_markers.get(g, {"p_value_ttest": 1, "effect_size_ttest": 0})
                b_stats.update(t_stats)
                b_stats.update(
                    {
                        "tissue_ontology_term_id": tiss,
                        "organism_ontology_term_id": organism,
                        "cell_type_ontology_term_id": ct,
                        "gene_ontology_term_id": g,
                    }
                )
                markers.append(b_stats)

            df = pd.DataFrame.from_records(markers)
            if df.shape[0] > 0:
                df = df.astype(
                    {
                        "p_value_binomtest": "float32",
                        "effect_size_binomtest": "float32",
                        "p_value_ttest": "float32",
                        "effect_size_ttest": "float32",
                    }
                )
                tiledb.from_pandas(uri, df, mode="append")

    logger.debug("Cube created, start consolidation")
    tiledb.consolidate(uri)

    logger.debug("Cube consolidated, start vacuumming")
    tiledb.vacuum(uri)

    logger.info(f"Marker genes cube created and stored at {uri}")
