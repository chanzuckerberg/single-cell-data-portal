import os

ONTOLOGY_TREE_FILENAME = "ontology_graph.json"
CELL_GUIDE_METADATA_FILENAME = "celltype_metadata.json"
CELL_GUIDE_TISSUE_METADATA_FILENAME = "tissue_metadata.json"
ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME = "cell_type_ontology_tree_state"
ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME = "tissue_ontology_tree_state"
CANONICAL_MARKER_GENES_FOLDERNAME = "canonical_marker_genes"
COMPUTATIONAL_MARKER_GENES_FOLDERNAME = "computational_marker_genes"
SOURCE_COLLECTIONS_FOLDERNAME = "source_collections"
GPT_OUTPUT_DIRECTORY_FOLDERNAME = "gpt_descriptions"

UBERON_BASIC_PERMANENT_URL_PRONTO = "http://purl.obolibrary.org/obo/uberon.obo"

ASCTB_MASTER_SHEET_URL = "https://purl.org/ccf/releases/2.2.1/ccf-asctb-all.json"

ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = "ensembl_gene_ids_to_descriptions.tsv.gz"

HOMO_SAPIENS_ORGANISM_ONTOLOGY_TERM_ID = "NCBITaxon:9606"

# If CELLGUIDE_PIPELINE_NUM_CPUS is not set, use 24 CPUs by default
# If the number of CPUs on the machine is less than 24, use the number of CPUs on the machine
# 24 CPUs was chosen to balance memory usage and speed on a c6i.32xlarge EC2 machine
# In trial runs, the memory usage did not exceed 50% of the available memory, which provides
# ample buffer.
CELLGUIDE_PIPELINE_NUM_CPUS = min(os.cpu_count(), os.getenv("CELLGUIDE_PIPELINE_NUM_CPUS", 12))

CELL_GUIDE_DATA_BUCKET_PATH_PREFIX = "s3://cellguide-data-public-"
