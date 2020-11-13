import collections


class CorporaConstants(object):
    # Namedtuples used to encapsulate information about a piece of metadata for verification purposes
    TypedMetadata = collections.namedtuple(
        "TypedMetadata", "field_name required_type valid_alternative", defaults=[None]
    )
    Ontology = collections.namedtuple("Ontology", "ontology_name s3_uri")

    # Constants related to type checking ontology fields
    EFO_ONTOLOGY = Ontology(ontology_name="efo", s3_uri="corpora-data-portal-ontologies/efo.txt")
    CELL_ONTOLOGY = Ontology(ontology_name="cl", s3_uri="corpora-data-portal-ontologies/cl.txt")
    HANCESTRO_ONTOLOGY = Ontology(ontology_name="hancestro", s3_uri="corpora-data-portal-ontologies/hancestro.txt")
    HSAPDV_ONTOLOGY = Ontology(ontology_name="hsapdv", s3_uri="corpora-data-portal-ontologies/hsapdv.txt")
    MONDO_ONTOLOGY = Ontology(ontology_name="mondo", s3_uri="corpora-data-portal-ontologies/mondo.txt")
    NCBI_TAXON_ONTOLOGY = Ontology(ontology_name="ncbi_taxon", s3_uri="corpora-data-portal-ontologies/ncbitaxon.txt")
    UBERON_ONTOLOGY = Ontology(ontology_name="uberon", s3_uri="corpora-data-portal-ontologies/uberon.txt")

    CORPORA_ONTOLOGIES = [
        EFO_ONTOLOGY,
        CELL_ONTOLOGY,
        HANCESTRO_ONTOLOGY,
        HSAPDV_ONTOLOGY,
        MONDO_ONTOLOGY,
        NCBI_TAXON_ONTOLOGY,
        UBERON_ONTOLOGY,
    ]

    # Constants related to type checking enum fields
    SEX_VALID_ENUM_VALUES = ["male", "female", "na", "unknown"]

    # Constants related to processing layers
    LAYER_DESCRIPTIONS = "layer_descriptions"
    X_DATA_LAYER_NAME = "X"
    RAW_DATA_LAYER_NAME = "raw.X"

    # Constants related to the metadata requirements of submitting a dataset to Corpora. Note that these metadata
    # fields are CASE-SENSITIVE.
    REQUIRED_OBSERVATION_METADATA_FIELDS = [
        TypedMetadata(field_name="tissue", required_type=str),
        TypedMetadata(field_name="assay", required_type=str),
        TypedMetadata(field_name="disease", required_type=str),
        TypedMetadata(field_name="cell_type", required_type=str),
        TypedMetadata(field_name="sex", required_type=SEX_VALID_ENUM_VALUES),
        TypedMetadata(field_name="ethnicity", required_type=str),
        TypedMetadata(field_name="developmental_stage", required_type=str),
    ]
    REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS = [
        TypedMetadata(field_name="tissue_ontology", required_type=UBERON_ONTOLOGY),
        TypedMetadata(field_name="assay_ontology", required_type=EFO_ONTOLOGY),
        TypedMetadata(field_name="disease_ontology", required_type=MONDO_ONTOLOGY, valid_alternative="PATO:0000461"),
        TypedMetadata(field_name="cell_type_ontology", required_type=CELL_ONTOLOGY, valid_alternative=""),
        TypedMetadata(field_name="ethnicity_ontology", required_type=HANCESTRO_ONTOLOGY, valid_alternative="na"),
        TypedMetadata(
            field_name="developmental_stage_ontology", required_type=HSAPDV_ONTOLOGY, valid_alternative="EFO:0000399"
        ),
    ]
    REQUIRED_DATASET_METADATA_FIELDS = [
        TypedMetadata(field_name="organism", required_type=str),
        TypedMetadata(field_name="organism_ontology", required_type=NCBI_TAXON_ONTOLOGY),
        TypedMetadata(field_name=LAYER_DESCRIPTIONS, required_type=dict),
    ]
    REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS = [
        TypedMetadata(field_name="title", required_type=str),
    ]

    REQUIRED_DATASET_PRESENTATION_HINTS_METADATA_FIELDS = [
        TypedMetadata(field_name="color_map", required_type=dict),
        TypedMetadata(field_name="tags", required_type=list),
        TypedMetadata(field_name="default_field", required_type=str),
        TypedMetadata(field_name="default_embedding", required_type=str),
    ]
    OPTIONAL_COLLECTION_LEVEL_METADATA_FIELDS = [
        TypedMetadata(field_name="collection_name", required_type=str),
        TypedMetadata(field_name="collection_description", required_type=str),
        TypedMetadata(field_name="collection_protocol_links", required_type=list),
        TypedMetadata(field_name="collection_raw_data_links", required_type=list),
        TypedMetadata(field_name="collection_other_links", required_type=list),
    ]

    # Constants related to the dataset format of which a submission to Corpora must be.
    LOOM_FILE_TYPE = ".loom"
    SEURAT_FILE_TYPE = ".rds"
    H5AD_FILE_TYPE = ".h5ad"
