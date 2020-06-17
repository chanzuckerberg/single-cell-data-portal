class CorporaConstants(object):
    # Constants related to processing layers
    LAYERS_DESCRIPTIONS = "layer_descriptions"
    X_DATA_LAYER_NAME = "x"
    RAW_DATA_LAYER_NAME = "raw.x"

    # Constants related to the metadata requirements of submitting a dataset to Corpora. Note that these metadata
    # fields are CASE-SENSITIVE.
    REQUIRED_OBSERVATION_METADATA_FIELDS = ["tissue", "assay", "disease", "cell_type", "sex", "ethnicity",
                                            "developmental_stage"]
    REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS = [
        "tissue_ontology",
        "assay_ontology",
        "disease_ontology",
        "cell_type_ontology",
        "ethnicity_ontology",
        "developmental_stage_ontology",
    ]
    REQUIRED_DATASET_METADATA_FIELDS = ["organism", "organism_ontology", LAYERS_DESCRIPTIONS]
    REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS = ["title", "contributors", "preprint_doi", "publication_doi"]
    REQUIRED_DATASET_PRESENTATION_HINTS_METADATA_FIELDS = ["color_map", "tags", "default_field", "default_embedding"]
    OPTIONAL_PROJECT_LEVEL_METADATA_FIELDS = [
        "project_name",
        "project_description",
        "project_protocol_links",
        "project_raw_data_links",
        "project_other_links",
    ]

    # Constants related to the dataset format of which a submission to Corpora must be.
    LOOM_FILE_TYPE = ".loom"
    SEURAT_FILE_TYPE = ".rds"
    H5AD_FILE_TYPE = ".h5ad"
