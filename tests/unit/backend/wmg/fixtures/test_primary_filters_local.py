test_organism_terms = [
    {"organism_ontology_term_id_0": "Mus musculus"},
    {"organism_ontology_term_id_1": "Callithrix jacchus"},
    {"organism_ontology_term_id_2": "Homo sapiens"}
]

test_tissue_terms = {
    "organism_ontology_term_id_0": [
        {"tissue_ontology_term_id_0": "primary motor cortex"},
        {"tissue_ontology_term_id_1": "urethra"},
        {"tissue_ontology_term_id_2": "large intestine"},
    ],
    "organism_ontology_term_id_1": [
        {"tissue_ontology_term_id_0": "primary motor cortex"},
        {"tissue_ontology_term_id_1": "urethra"},
        {"tissue_ontology_term_id_2": "large intestine"},
    ],
    "organism_ontology_term_id_2": [
        {"tissue_ontology_term_id_0": "primary motor cortex"},
        {"tissue_ontology_term_id_1": "urethra"},
        {"tissue_ontology_term_id_2": "large intestine"},
    ],
}
test_gene_terms = {
    "organism_ontology_term_id_0": [
        {"gene_ontology_term_id_0": "TSPAN6"},
        {"gene_ontology_term_id_1": "TNMD"},
        {"gene_ontology_term_id_2": "SCYL3"},
    ],
    "organism_ontology_term_id_1": [
        {"gene_ontology_term_id_0": "TSPAN6"},
        {"gene_ontology_term_id_1": "TNMD"},
        {"gene_ontology_term_id_2": "SCYL3"},
    ],
    "organism_ontology_term_id_2": [
        {"gene_ontology_term_id_0": "TSPAN6"},
        {"gene_ontology_term_id_1": "TNMD"},
        {"gene_ontology_term_id_2": "SCYL3"},
    ],        
}

test_snapshot_id = "test-snapshot-id"


def build_precomputed_primary_filters_local(
    snapshot_id=test_snapshot_id,
    organism_terms=test_organism_terms,
    tissue_terms=test_tissue_terms,
    gene_terms=test_gene_terms,
):

    return {
        "snapshot_id": snapshot_id,
        "organism_terms": organism_terms,
        "tissue_terms": tissue_terms,
        "gene_terms": gene_terms,
    }
