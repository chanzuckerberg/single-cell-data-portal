test_organism_terms = [
    {"NCBITaxon:10090": "Mus musculus"},
    {"NCBITaxon:9483": "Callithrix jacchus"},
    {"NCBITaxon:9544": "Macaca mulatta"},
    {"NCBITaxon:9606": "Homo sapiens"},
]

test_tissue_terms = {
    "NCBITaxon:10090": [
        {"UBERON:0001384": "primary motor cortex"},
        {"UBERON:0000057": "urethra"},
        {"UBERON:0000059": "large intestine"},
        {"UBERON:0000948": "heart"},
        {"UBERON:0001013": "adipose tissue"},
    ],
    "NCBITaxon:9483": [{"UBERON:0001384": "primary motor cortex"}],
    "NCBITaxon:9544": [{"UBERON:0001384": "primary motor cortex"}, {"UBERON:0002048": "lung"}],
    "NCBITaxon:9606": [
        {"UBERON:0000029": "lymph node"},
        {"UBERON:0000115": "lung epithelium"},
        {"UBERON:0000160": "intestine"},
        {"UBERON:0000178": "blood"},
        {"UBERON:0000362": "renal medulla"},
    ],
}
test_gene_terms = {
    "NCBITaxon:10090": [
        {"ENSG00000000003": "TSPAN6"},
        {"ENSG00000000005": "TNMD"},
        {"ENSG00000000457": "SCYL3"},
        {"ENSG00000000460": "C1orf112"},
        {"ENSG00000000938": "FGR"},
    ],
    "NCBITaxon:9483": [
        {"ENSG00000000003": "TSPAN6"},
        {"ENSG00000000005": "TNMD"},
        {"ENSG00000000419": "DPM1"},
        {"ENSG00000000457": "SCYL3"},
        {"ENSG00000000460": "C1orf112"},
    ],
    "NCBITaxon:9544": [
        {"ENSG00000000419": "DPM1"},
        {"ENSG00000000457": "SCYL3"},
        {"ENSG00000000460": "C1orf112"},
        {"ENSG00000000938": "FGR"},
        {"ENSG00000001036": "FUCA2"},
    ],
    "NCBITaxon:9606": [
        {"ENSG00000000003": "TSPAN6"},
        {"ENSG00000000005": "TNMD"},
        {"ENSG00000000419": "DPM1"},
        {"ENSG00000000457": "SCYL3"},
        {"ENSG00000000460": "C1orf112"},
    ],
}

test_snapshot_id = "test-snapshot-id"


def build_precomputed_primary_filters(
    snapshot_id=test_snapshot_id,
    organism_terms=test_organism_terms,
    tissue_terms=test_tissue_terms,
    gene_terms=test_gene_terms,
):
    return {
        snapshot_id: snapshot_id,
        organism_terms: organism_terms,
        tissue_terms: tissue_terms,
        gene_terms: gene_terms,
    }
