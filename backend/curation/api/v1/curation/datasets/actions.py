from flask import jsonify, make_response

from backend.curation.api.v1.curation.collections.common import reshape_dataset_for_curation_api
from backend.portal.api.providers import get_business_logic


def get():
    """
    Datasets index endpoint to retrieve full metadata for all public Datasets.
    """
    datasets = get_business_logic().get_all_published_datasets()

    return make_response(
        jsonify([reshape_dataset_for_curation_api(d, is_published=True, use_canonical_id=True) for d in datasets]), 200
    )
