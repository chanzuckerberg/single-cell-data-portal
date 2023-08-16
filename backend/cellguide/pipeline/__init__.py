import sys
import time

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline


def run_cellguide_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Run ontology tree pipeline
    ontology_tree = run_ontology_tree_pipeline(output_directory)

    run_canonical_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)


if __name__ == "__main__":
    run_cellguide_pipeline()
    sys.exit()
