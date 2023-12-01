"""
This method enriches the response of a Dataset object with ancestor mappings.
These responses are specific to the current API layer
"""

from collections import OrderedDict

from backend.common.feature_flag import FeatureFlagService, FeatureFlagValues


def enrich_dataset_with_ancestors(dataset, key, ontology_mapping):
    """
    Tag dataset with ancestors for all values of the given key, if any.
    """
    if key not in dataset:
        return

    terms = [e["ontology_term_id"] for e in dataset[key]]

    is_schema_4 = FeatureFlagService.is_enabled(FeatureFlagValues.SCHEMA_4)
    is_tissue = key == "tissue"
    if is_tissue and is_schema_4:
        # TODO remove is_schema_4 condition once Schema 4 is rolled out and
        # feature flag is removed (#6266). "tissue" must include "tissue_type"
        # when generating ancestors; "cell_type" and "development_stage" do not.
        terms = [generate_tagged_tissue_ontology_id(e) for e in dataset[key]]
    else:
        terms = [e["ontology_term_id"] for e in dataset[key]]

    if not terms:
        return

    ancestors = [ontology_mapping.get(term) for term in terms]
    flattened_ancestors = [item for sublist in ancestors if sublist for item in sublist]
    unique_ancestors = list(OrderedDict.fromkeys(flattened_ancestors))
    if unique_ancestors:
        dataset[f"{key}_ancestors"] = unique_ancestors


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
