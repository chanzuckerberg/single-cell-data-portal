import os

ONTOLOGY_TREE_FILENAME = "ontology_graph.json"
CELL_GUIDE_METADATA_FILENAME = "celltype_metadata.json"
CELL_GUIDE_VALID_EXPLORER_CXGS_FILENAME = "valid_explorer_cxgs.json"
CELL_GUIDE_CELL_TYPE_MAPPING_FILENAME = "celltype_to_tissue_mapping.json"
ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME = "ontology_tree"
CELL_GUIDE_TISSUE_METADATA_FILENAME = "tissue_metadata.json"
ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME = "cell_type_ontology_tree_state"
ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME = "tissue_ontology_tree_state"
CANONICAL_MARKER_GENES_FOLDERNAME = "canonical_marker_genes"
COMPUTATIONAL_MARKER_GENES_FOLDERNAME = "computational_marker_genes"
SOURCE_COLLECTIONS_FOLDERNAME = "source_collections"
GPT_OUTPUT_DIRECTORY_FOLDERNAME = "gpt_descriptions"
GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME = "gpt_seo_descriptions"
MARKER_GENE_PRESENCE_FILENAME = "marker_gene_presence.json.gz"


ASCTB_MASTER_SHEET_URL = "https://ccf-ontology.hubmapconsortium.org/v2.3.0/ccf-asctb-all.json"

ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = "ensembl_gene_ids_to_descriptions.tsv.gz"

HOMO_SAPIENS_ORGANISM_ONTOLOGY_TERM_ID = "NCBITaxon:9606"

# If CELLGUIDE_PIPELINE_NUM_CPUS is not set, use 12 CPUs by default
# If the number of CPUs on the machine is less than 12, use the number of CPUs on the machine
# TODO: Tune this number. But note - speed is not that important here and we've run into issues
# where the number of CPUs that was previously set (24) became too high and resulted in OOM isues
# as the data corpus grew.
CELLGUIDE_PIPELINE_NUM_CPUS = min(os.cpu_count(), os.getenv("CELLGUIDE_PIPELINE_NUM_CPUS", 12))

CELL_GUIDE_DATA_BUCKET_PATH_PREFIX = "s3://cellguide-data-public-"

CELL_GUIDE_PINNED_SCHEMA_VERSION = "5.0.0"
