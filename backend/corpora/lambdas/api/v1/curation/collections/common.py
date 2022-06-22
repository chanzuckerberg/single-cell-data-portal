from ...authorization import is_user_owner_or_allowed
from ......common.corpora_config import CorporaConfig
from ......common.corpora_orm import CollectionVisibility


def reshape_for_curation_api_and_is_allowed(collection, token_info, allow_access=False):
    owner = collection["owner"]
    if is_user_owner_or_allowed(token_info, owner):
        collection["access_type"] = "WRITE"
    elif not allow_access and collection["visibility"] == CollectionVisibility.PRIVATE:
        # User neither provided the uuid for access nor are they authorized by their access token
        return False
    elif token_info:
        # Access token was provided but user is not authorized
        collection["access_type"] = "READ"
    else:
        # No access token was provided
        collection["access_type"] = None

    del collection["owner"]  # Don't actually want to return 'owner' in response
    collection["collection_url"] = f"{CorporaConfig().collections_base_url}/collections/{collection['id']}"

    if "datasets" in collection:
        for dataset in collection["datasets"]:
            if "artifacts" in dataset:
                dataset["dataset_assets"] = dataset.pop("artifacts")
            if "processing_status" in dataset:
                dataset["processing_status"] = dataset["processing_status"]["processing_status"]

    return True
