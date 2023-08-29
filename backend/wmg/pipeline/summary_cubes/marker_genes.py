import contextlib
import gc
import logging

import numpy as np
import pandas as pd
import tiledb

from backend.common.utils.exceptions import MarkerGeneCalculationException
from backend.wmg.data.schemas.marker_gene_cube_schema import marker_genes_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, MARKER_GENES_CUBE_NAME
from backend.wmg.data.utils import create_empty_cube, log_func_runtime
from backend.wmg.pipeline.summary_cubes.calculate_markers import get_markers

logger: logging.Logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@contextlib.contextmanager
def extract_tissue_celltype_organism(corpus_path: str) -> pd.DataFrame:
    # extract cell counts grouped by tissue, cell type, and organism
    with tiledb.open(f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}") as array:
        yield (
            array.query(
                attrs=["cell_type_ontology_term_id"],
                dims=["organism_ontology_term_id", "tissue_ontology_term_id"],
            )
            .df[:]
            .groupby(
                [
                    "tissue_ontology_term_id",
                    "cell_type_ontology_term_id",
                    "organism_ontology_term_id",
                ]
            )
            .first()
        )


@log_func_runtime
def create_marker_genes_cube(corpus_path: str):
    """
    Create marker genes cube and write to disk
    """
    with extract_tissue_celltype_organism(corpus_path) as tco:
        uniq_tissues = tco.index.levels[0]
        tissues = np.array(list(tco.index.get_level_values(0)))
        cell_types = np.array(list(tco.index.get_level_values(1)))
        organisms = np.array(list(tco.index.get_level_values(2)))

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
            try:
                t_markers = get_markers(
                    target, context, corpus=corpus_path, test="ttest", percentile=0.05, n_markers=None
                )
            except MarkerGeneCalculationException as e:
                # exception handling here so pipeline doesn't fail if no cells match query criteria
                logger.info("Error finding markers for tissue: %s, cell type: %s, organism: %s", tiss, ct, organism)
                logger.info(e)
                continue
            gc.collect()

            all_marker_genes = set(t_markers.keys())
            markers = []
            for g in all_marker_genes:
                t_stats = t_markers.get(g, {"p_value_ttest": np.nan, "effect_size_ttest": np.nan})
                t_stats.update(
                    {
                        "tissue_ontology_term_id": tiss,
                        "organism_ontology_term_id": organism,
                        "cell_type_ontology_term_id": ct,
                        "gene_ontology_term_id": g,
                    }
                )
                markers.append(t_stats)

            df = pd.DataFrame.from_records(markers)
            if df.shape[0] > 0:
                df = df.astype(
                    {
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
