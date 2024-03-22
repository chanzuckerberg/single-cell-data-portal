"""
This method enriches the response of a Dataset object with ancestor mappings.
These responses are specific to the current API layer
"""

from cellxgene_ontology_guide.curated_ontology_term_lists import get_curated_ontology_term_list
from cellxgene_ontology_guide.entities import CuratedOntologyTermList
from cellxgene_ontology_guide.ontology_parser import OntologyParser

ONTOLOGY_PARSER = OntologyParser(schema_version="v5.0.0")  # TODO: this should be a constant
# TODO: clean-up where to set these constants
ACCEPTED_TISSUE_ANCESTORS = {
    term
    for term_list in ONTOLOGY_PARSER.map_term_descendants(
        get_curated_ontology_term_list(CuratedOntologyTermList.SYSTEM), include_self=True
    ).values()
    for term in term_list
}
ACCEPTED_CELL_TYPE_ANCESTORS = {
    term
    for term_list in ONTOLOGY_PARSER.map_term_descendants(
        get_curated_ontology_term_list(CuratedOntologyTermList.CELL_CLASS), include_self=True
    ).values()
    for term in term_list
}
ACCEPTED_UBERON_DEVELOPMENT_STAGE_ANCESTORS = {
    term
    for term_list in ONTOLOGY_PARSER.map_term_descendants(
        get_curated_ontology_term_list(CuratedOntologyTermList.UBERON_DEVELOPMENT_STAGE), include_self=True
    ).values()
    for term in term_list
}


def enrich_dataset_with_ancestors(dataset, key):
    """
    Tag dataset with ancestors for all values of the given key, if any.
    """
    if key not in dataset:
        return
    is_tissue = key == "tissue"
    unique_ancestors = set()
    for term in dataset[key]:
        ancestors = ONTOLOGY_PARSER.get_term_ancestors(term["ontology_term_id"], include_self=False)
        # TODO: clean-up this logic
        if is_tissue:
            ancestors = list(set(ancestors) & ACCEPTED_TISSUE_ANCESTORS)
        elif key == "cell_type":
            ancestors = list(set(ancestors) & ACCEPTED_CELL_TYPE_ANCESTORS)
        elif key == "development_stage" and "UBERON" in term["ontology_term_id"]:
            ancestors = list(set(ancestors) & ACCEPTED_UBERON_DEVELOPMENT_STAGE_ANCESTORS)
        # If the term is a tissue, tag itself with the tissue type in the ancestor list
        self_as_ancestor = generate_tagged_tissue_ontology_id(term) if is_tissue else term["ontology_term_id"]
        ancestors.append(self_as_ancestor)
        unique_ancestors.update(ancestors)
    if unique_ancestors:
        dataset[f"{key}_ancestors"] = list(unique_ancestors)


def generate_tagged_tissue_ontology_id(tissue):
    """
    Generate ontology ID tagged with tissue_type for the given tissue. For
    example, UBERON:1234567 (organoid).
    """
    tissue_id = tissue["ontology_term_id"]
    # Handle possible None for tissue_type (possible during migration): default
    # to "tissue".
    tissue_type = tissue["tissue_type"] or "tissue"
    if tissue_type == "tissue":
        return tissue_id
    return f"{tissue_id} ({tissue_type})"
