from backend.curation.api.v1.curation.collections.common import (
    validate_uuid_else_forbidden,
)
from backend.portal.api.providers import get_business_logic


def get(dataset_version_id: str):
    """
    Fetch Dataset metadata for a Dataset version.
    """
    business_logic = get_business_logic()
    validate_uuid_else_forbidden(dataset_version_id)

    business_logic.get_dataset_version(dataset_version_id)
    # get collection version ID
    pass
