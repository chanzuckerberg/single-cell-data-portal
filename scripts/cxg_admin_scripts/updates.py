import logging
import os
import sys
import time

import requests
from furl import furl

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

logging.basicConfig()
logger: logging.Logger = logging.getLogger(__name__)

auth0_apis = {
    "staging": "https://czi-cellxgene-dev.us.auth0.com",
    "dev": "https://czi-cellxgene-dev.us.auth0.com",
    "prod": "https://corpora-prod.auth0.com",
}


def add_trailing_slash_to_explorer_urls(ctx):
    """
    The explorer_url for datasets must end with a trailing slash to function
    properly. This script adds a trailing slash to a dataset's explorer_url
    if it already does not end with one.
    """
    pass


def get_collections_without_curator():
    pass


def get_owner_info_from_auth0(owner, access_token, deployment):
    auth0_api = auth0_apis[deployment]

    logger.info(f"Retrieving curator info for owner:{owner} from Auth0.")
    url_manager = furl(
        f"{auth0_api}/api/v2/users/{owner}",
        query_params=dict(fields=",".join(["given_name", "family_name", "name"]), include_fields="true"),
    )
    response = requests.get(url_manager.url, headers={"Authorization": f"Bearer {access_token}"})
    response.raise_for_status()
    if response.headers["X-RateLimit-Remaining"] == 0:
        wait = time.time() - response.headers["X-RateLimit-Reset"]
        time.sleep(wait)
    body = response.json()
    name = f"{body.get('given_name', '')} {body.get('family_name', '')}".strip() or body.get("name", "")
    return name


def update_database_curator_name(owner_id, owner_name):
    pass


def update_curator_names(ctx, access_token):
    """Add the curator name to all collection based on the owner of the collection.

    ACCESS_TOKEN: Retrieved from Auth0 console or generated using the Client ID and Client Secret.
    The application must be authorized to access to the Auth0 Management API with the following permissions read:users
    read:user_idp_tokens.
    """
    pass


def add_publisher_metadata(ctx):
    """Add publisher metadata to the current records"""
    pass


def refresh_preprint_doi(ctx):
    """Add publisher metadata to the current records"""
    pass


def update_collection_owner(ctx, collection_id, new_owner):
    """Update the owner of a cellxgene collection.
    To run:
    ./scripts/cxg_admin.py --deployment prod update-collection-owner "$COLLECTION_ID $NEW_OWNER_ID
    """
    pass


def transfer_collections(ctx, curr_owner, new_owner):
    """Transfer all collections owned by the curr_owner to the new_owner.
    Retrieve user ids from auth0 before running or ping an engineer on the team to check the id of the owner in the
    database
    To run:
    ./scripts/cxg_admin.py --deployment prod transfer-collections $CURR_OWNER_ID $NEW_OWNER_ID
    """
    pass


def strip_all_collection_fields(ctx):
    """
    Strip all the `collection` string fields, so whitespace at the beginning and the end are removed.
    """
    pass
