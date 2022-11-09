"""
This method enriches the response of a Dataset object with ancestor mappings.
These responses are specific to the current API layer
"""

from collections import OrderedDict


def enrich_dataset_with_ancestors(dataset, key, ontology_mapping):
    """
    Tag dataset with ancestors for all values of the given key, if any.
    """
    if key not in dataset:
        return

    terms = [e["ontology_term_id"] for e in dataset[key]]

    if not terms:
        return

    ancestors = [ontology_mapping.get(term) for term in terms]
    flattened_ancestors = [item for sublist in ancestors if sublist for item in sublist]
    unique_ancestors = list(OrderedDict.fromkeys(flattened_ancestors))
    if unique_ancestors:
        dataset[f"{key}_ancestors"] = unique_ancestors