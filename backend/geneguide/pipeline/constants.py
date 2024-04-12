import os

ONTOLOGY_TREE_FILENAME = "ontology_graph.json"
GENEGUIDE_METADATA_FILENAME = "go_term_metadata.json"
GENEGUIDE_VALID_EXPLORER_CXGS_FILENAME = "valid_explorer_cxgs.json"
ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME = "ontology_tree"
ONTOLOGY_TREE_STATE_PER_GOTERM_FOLDERNAME = "go_term_tree_state"


CANONICAL_MARKER_GENES_FOLDERNAME = "canonical_marker_genes"
COMPUTATIONAL_MARKER_GENES_FOLDERNAME = "computational_marker_genes"
SOURCE_COLLECTIONS_FOLDERNAME = "source_collections"
GPT_OUTPUT_DIRECTORY_FOLDERNAME = "gpt_descriptions"
GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME = "gpt_seo_descriptions"
MARKER_GENE_PRESENCE_FILENAME = "marker_gene_presence.json.gz"

UBERON_BASIC_PERMANENT_URL_PRONTO = "http://purl.obolibrary.org/obo/uberon.obo"

ASCTB_MASTER_SHEET_URL = "https://ccf-ontology.hubmapconsortium.org/v2.3.0/ccf-asctb-all.json"

ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = "ensembl_gene_ids_to_descriptions.tsv.gz"

HOMO_SAPIENS_ORGANISM_ONTOLOGY_TERM_ID = "NCBITaxon:9606"

# If GENEGUIDE_PIPELINE_NUM_CPUS is not set, use 24 CPUs by default
# If the number of CPUs on the machine is less than 24, use the number of CPUs on the machine
# 24 CPUs was chosen to balance memory usage and speed on a c6i.32xlarge EC2 machine
# In trial runs, the memory usage did not exceed 50% of the available memory, which provides
# ample buffer.
GENEGUIDE_PIPELINE_NUM_CPUS = min(os.cpu_count(), os.getenv("GENEGUIDE_PIPELINE_NUM_CPUS", 24))

GENEGUIDE_DATA_BUCKET_PATH_PREFIX = "s3://cellguide-data-public-"

GO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
