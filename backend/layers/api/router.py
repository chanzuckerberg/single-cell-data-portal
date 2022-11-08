# TODO: methods are temporarily defined here. Eventually, this should route calls to an appropriate class


from typing import Optional
# from backend.layers.api.portal_api import PortalApi

def portal_api():
    # TODO: put the real PortalApi here. Tests will mock this method to return a testable version of PortalAPI.
    # This is a poor man dependency injection, basically
    return NotImplemented

def get_collections_list(from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
    return portal_api().get_collections_list(from_date, to_date, token_info)

def get_collection_details(collection_id: str, token_info: dict):
    return portal_api().get_collection_details(collection_id, token_info)

def get_collections_index():
    return portal_api().get_collection_index()

def post_collection_revision(collection_id: str, token_info: dict):
    return portal_api().post_collection_revision(collection_id, token_info)

def create_collection(body: dict, user: str):
    return portal_api().create_collection(body, user)

def delete_collection(collection_id: str, token_info: dict):
    return portal_api().delete_collection(collection_id, token_info)

def update_collection(collection_id: str, body: dict, token_info: dict):
    return portal_api().update_collection(collection_id, body, token_info)


def post(collection_id: str, body: object, token_info: dict): # publish
    pass

def link(collection_id: str, body: dict, token_info: dict):
    pass

def relink(collection_id: str, body: dict, token_info: dict):
    pass

def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None):
    pass



def post_dataset_asset(dataset_id: str, asset_id: str):
    pass

def get_dataset_assets(dataset_id: str):
    pass

def get_status(dataset_id: str, token_info: dict):
    pass

def get_datasets_index():
    pass

def delete_dataset(dataset_id: str, token_info: dict):
    pass

def get_dataset_identifiers(url: str):
    pass

def post_dataset_gene_sets(dataset_id: str, body: object, token_info: dict):
    pass
