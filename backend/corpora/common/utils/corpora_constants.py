class CorporaConstants(object):
    # Constants related to processing layers
    LAYERS_DESCRIPTIONS = "LAYERS_DESCRIPTIONS"
    X_DATA_LAYER_NAME = "X"
    RAW_DATA_LAYER_NAME = "RAW.X"

    # Constants related to the metadata requirements of submitting a dataset to Corpora.
    REQUIRED_OBSERVATION_METADATA_FIELDS = [
        "TISSUE",
        "ASSAY",
        "DISEASE",
        "CELL_TYPE",
        "SEX",
        "ETHNICITY",
        "DEVELOPMENTAL_STAGE",
    ]
    REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS = [
        "TISSUE_ONTOLOGY",
        "ASSAY_ONTOLOGY",
        "DISEASE_ONTOLOGY",
        "CELL_TYPE_ONTOLOGY",
        "ETHNICITY_ONTOLOGY",
        "DEVELOPMENTAL_STAGE_ONTOLOGY",
    ]
    REQUIRED_DATASET_METADATA_FIELDS = ["ORGANISM", "ORGANISM_ONTOLOGY", LAYERS_DESCRIPTIONS]
    REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS = ["TITLE", "CONTRIBUTORS", "PREPRINT_DOI", "PUBLICATION_DOI"]
    REQUIRED_DATASET_PRESENTATION_HINTS_METADATA_FIELDS = ["COLOR_MAP", "TAGS", "DEFAULT_FIELD", "DEFAULT_EMBEDDING"]
    OPTIONAL_PROJECT_LEVEL_METADATA_FIELDS = [
        "PROJECT_NAME",
        "PROJECT_DESCRIPTION",
        "PROJECT_PROTOCOL_LINKS",
        "PROJECT_RAW_DATA_LINKS",
        "PROJECT_OTHER_LINKS",
    ]

    # Constants related to the dataset format of which a submission to Corpora must be.
    LOOM_FILE_TYPE = ".loom"
    SEURAT_FILE_TYPE = ".rds"
    H5AD_FILE_TYPE = ".h5ad"
