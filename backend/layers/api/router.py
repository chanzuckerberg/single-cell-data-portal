# TODO: methods are temporarily defined here. Eventually, this should route calls to an appropriate class


from typing import Optional

from backend.layers.api.portal_api import PortalApi
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.crossref_provider import CrossrefProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProvider

api = None
_business_logic = None


def get_business_logic():
    global _business_logic
    if not _business_logic:
        _business_logic = BusinessLogic(
            DatabaseProvider(), CrossrefProvider(), StepFunctionProvider(), S3Provider(), UriProvider()
        )
    return _business_logic


# TODO: put the real PortalApi here. Tests will mock this method to return a testable version of PortalAPI.
# This is a poor man dependency injection, basically
def portal_api():
    global api
    if api is None:
        database_provider = DatabaseProvider()
        crossref_provider = CrossrefProvider()
        step_function_provider = StepFunctionProvider()
        s3_provider = S3Provider()
        uri_provider = UriProvider()

        business_logic = BusinessLogic(
            database_provider, crossref_provider, step_function_provider, s3_provider, uri_provider
        )
        portal_api = PortalApi(business_logic)
        api = portal_api
    return api


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


def publish_post(collection_id: str, body: object, token_info: dict):
    return portal_api().publish_post(collection_id, body, token_info)


def upload_link(collection_id: str, body: dict, token_info: dict):
    return portal_api().upload_link(collection_id, body, token_info)


def upload_relink(collection_id: str, body: dict, token_info: dict):
    return portal_api().upload_relink(collection_id, body, token_info)


def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None):
    return portal_api().upload_from_link(collection_id, token_info, url, dataset_id)


def post_dataset_asset(dataset_id: str, asset_id: str):
    return portal_api().post_dataset_asset(dataset_id, asset_id)


def get_dataset_assets(dataset_id: str):
    return portal_api().get_dataset_assets(dataset_id)


def get_status(dataset_id: str, token_info: dict):
    return portal_api().get_status(dataset_id, token_info)


def get_datasets_index():
    return portal_api().get_datasets_index()


def delete_dataset(dataset_id: str, token_info: dict):
    return portal_api().delete_dataset(dataset_id, token_info)


def get_dataset_identifiers(url: str):
    return portal_api().get_dataset_identifiers(url)
