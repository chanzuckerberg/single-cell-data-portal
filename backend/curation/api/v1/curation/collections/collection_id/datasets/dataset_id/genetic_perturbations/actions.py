from flask import jsonify, make_response

from backend.curation.api.v1.curation.collections.common import _get_collection_and_dataset
from backend.portal.api.providers import get_business_logic


def get(collection_id: str, dataset_id: str):
    _, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)
    gp_metadata = get_business_logic().get_dataset_genetic_perturbations(dataset_version.version_id)
    return make_response(jsonify({"genetic_perturbations": gp_metadata.to_dict() if gp_metadata else None}), 200)
