class CorporaConstants(object):
    # Constants related to the metadata requirements of submitting a dataset to Corpora. Note that these metadata
    # fields are CASE-SENSITIVE.
    REQUIRED_OBSERVATION_METADATA_FIELDS = ["tissue", "assay", "disease", "cell_type", "sex", "ethnicity"]
    REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS = [
        "tissue_ontology",
        "assay_ontology",
        "disease_ontology",
        "cell_type_ontology",
        "ethnicity_ontology",
    ]
    REQUIRED_DATASET_METADATA_FIELDS = ["organism", "organism_ontology", "layer_descriptions"]
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
