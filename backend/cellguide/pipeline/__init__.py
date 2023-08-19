import shutil
import time

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline
from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline
from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline


def run_cellguide_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Run ontology tree pipeline
    ontology_tree = run_ontology_tree_pipeline(output_directory)

    # Generate cell guide cards, synonyms, and descriptions
    run_metadata_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

    # Generate canonical marker genes from ASCT-B (HUBMAP)
    run_canonical_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

    # Generate source data for each cell type
    run_source_collections_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

    # Generate computational marker genes from the CZI corpus
    run_computational_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

    # zip up the results
    shutil.make_archive(output_directory, "zip", output_directory)
