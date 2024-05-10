import os

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

CELL_GUIDE_PINNED_SCHEMA_VERSION = "5.1.0"
