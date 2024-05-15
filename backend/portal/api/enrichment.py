"""
This method enriches the response of a Dataset object with ancestor mappings.
These responses are specific to the current API layer
"""

from typing import Dict, Set

from cellxgene_ontology_guide.curated_ontology_term_lists import get_curated_ontology_term_list
from cellxgene_ontology_guide.entities import CuratedOntologyTermList
from cellxgene_ontology_guide.ontology_parser import OntologyParser

from backend.layers.common.entities import DatasetVersion

ONTOLOGY_PARSER = OntologyParser()

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


def enrich_dataset_with_ancestors(dataset: DatasetVersion, key: str, corpus_term_set: Set[str]) -> None:
    """
    Tag dataset with ancestors for all values of the given key, if any.
    """
    if key not in dataset:
        return
    is_tissue = key == "tissue"
    is_cell_type = key == "cell_type"
    is_development_stage = key == "development_stage"
    unique_ancestors = set()
    accepted_ancestors = set()
    accepted_ancestors.update(corpus_term_set)
    if is_tissue:
        accepted_ancestors.update(ACCEPTED_TISSUE_ANCESTORS)
    elif is_cell_type:
        accepted_ancestors.update(ACCEPTED_CELL_TYPE_ANCESTORS)
    elif is_development_stage:
        accepted_ancestors.update(ACCEPTED_UBERON_DEVELOPMENT_STAGE_ANCESTORS)
    for term in dataset[key]:
        if term == "unknown":
            continue
        ancestors = ONTOLOGY_PARSER.get_term_ancestors(term["ontology_term_id"], include_self=False)
        # no filtering of ancestors for non-UBERON development stage terms
        # otherwise, filter ancestors to only include ancestors that are in the accepted_ancestors set or from the
        # corpus term set
        if is_tissue or is_cell_type or (is_development_stage and "UBERON" in term["ontology_term_id"]):
            ancestors = list(set(ancestors) & accepted_ancestors)
        # If the term is a tissue, tag itself with the tissue type in the ancestor list
        self_as_ancestor = generate_tagged_tissue_ontology_id(term) if is_tissue else term["ontology_term_id"]
        ancestors.append(self_as_ancestor)
        unique_ancestors.update(ancestors)
    if unique_ancestors:
        dataset[f"{key}_ancestors"] = list(unique_ancestors)


def generate_tagged_tissue_ontology_id(tissue: Dict[str, str]) -> str:
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
