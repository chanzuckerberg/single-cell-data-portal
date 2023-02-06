import logging
import os
import sys
import time
import traceback

import click
import requests
from furl import furl
from requests import HTTPError
from sqlalchemy import null

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa
from backend.common.corpora_orm import (
    CollectionLinkType,
    CollectionVisibility,
    DbCollection,
    DbCollectionLink,
    DbDataset,
)
from backend.common.entities.collection import Collection
from backend.common.providers import crossref_provider
from backend.common.utils.db_session import db_session_manager

logging.basicConfig()
logger = logging.getLogger(__name__)

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
    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will add a trailing slash to "
            "a dataset's explorer_url if it already does not end with one.",
            abort=True,
        )

        for record in session.query(DbDataset):
            dataset_id = record.id

            if record.explorer_url is None:
                logger.info(f"SKIPPING - Dataset does not have an explorer_url | dataset_id {dataset_id}")
                continue

            explorer_url = record.explorer_url.strip()
            if explorer_url[-1] == "/":
                logger.info(
                    f"SKIPPING - Dataset explorer_url already ends with trailing slash | dataset_id {dataset_id}"
                )
                continue

            logger.info(
                f"Adding trailing slash to dataset explorer_url | dataset_id {dataset_id} | "
                f"original explorer_url: {explorer_url}"
            )
            explorer_url_w_slash = explorer_url + "/"
            record.explorer_url = explorer_url_w_slash

        logger.info("----- Finished adding trailing slash to explorer_url for datasets! ----")


def get_collections_without_curator():
    logger.info("Gathering collections with no curator.")
    with db_session_manager() as session:
        _filter = (Collection.table.curator_name == null()) | (Collection.table.curator_name == "")
        _owners = [result.owner for result in session.query(Collection.table.owner).filter(_filter).distinct().all()]
    return _owners


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
    logger.info(f"Updating collections owned by {owner_id} with curator_name:{owner_name}.")
    with db_session_manager() as session:
        collections = session.query(DbCollection).filter(DbCollection.owner == owner_id).all()
        for collection in collections:
            collection.curator_name = owner_name
            if collection.visibility == CollectionVisibility.PUBLIC:
                collection.data_submission_policy_version = "2.0"


def update_curator_names(ctx, access_token):
    """Add the curator name to all collection based on the owner of the collection.

    ACCESS_TOKEN: Retrieved from Auth0 console or generated using the Client ID and Client Secret.
    The application must be authorized to access to the Auth0 Management API with the following permissions read:users
    read:user_idp_tokens.
    """
    deployment = ctx.obj["deployment"]

    unique_owners = get_collections_without_curator()
    owners = dict()
    bad_owners = []
    for owner in unique_owners:
        try:
            owner_name = get_owner_info_from_auth0(owner, access_token, deployment)
            if owner_name:
                owners[owner] = owner_name
        except HTTPError:
            bad_owners.append(owner)
            logger.exception(f"Failed to fetch Auth0 info for owner:{owner}")
    for owner_id, owner_name in owners.items():
        update_database_curator_name(owner_id, owner_name)


def add_publisher_metadata(ctx):
    """Add publisher metadata to the current records"""

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will populate publisher_metadata for all "
            "datasets. This will also do N calls to Crossref.",
            abort=True,
        )

        provider = crossref_provider.CrossrefProvider()

        for record in session.query(DbCollection):
            collection = Collection.get_collection(session, record.id, record.visibility)
            if not collection:
                continue
            collection_id = record.id
            dois = [link.link_url for link in collection.links if link.link_type == CollectionLinkType.DOI]
            normalized_doi = collection.get_doi()

            if record.publisher_metadata:
                print("Already has metadata, skipping...")
                continue

            if normalized_doi:
                if "/" not in normalized_doi:
                    continue

                try:
                    metadata = provider.fetch_metadata(normalized_doi)
                    record.publisher_metadata = metadata
                    print(collection_id, dois, normalized_doi)
                except Exception as e:
                    print(record.id, normalized_doi, record.name, e)
                    print(traceback.format_exc())
                    record.publisher_metadata = {"error": str(e.__cause__)}

                time.sleep(2)


def refresh_preprint_doi(ctx):
    """Add publisher metadata to the current records"""

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will populate publisher_metadata for all "
            "datasets. This will also do N calls to Crossref.",
            abort=True,
        )

        provider = crossref_provider.CrossrefProvider()

        for record in session.query(DbCollection):
            collection = Collection.get_collection(session, record.id, record.visibility)
            if not collection:
                continue
            collection_id = record.id
            normalized_doi = collection.get_doi()

            if record.publisher_metadata and record.publisher_metadata.get("is_preprint"):
                print(f"{normalized_doi}: is_preprint")

                try:
                    published_doi = provider.fetch_preprint_published_doi(normalized_doi)
                    print(collection_id, normalized_doi, published_doi)
                    doi_record = session.query(DbCollectionLink).filter(
                        DbCollectionLink.collection_id == record.id,
                        DbCollectionLink.collection_visibility == record.visibility,
                        DbCollectionLink.link_type == CollectionLinkType.DOI,
                    )
                    doi_record.link_url = published_doi
                except Exception as e:
                    print(e)


def update_collection_owner(ctx, collection_id, new_owner):
    """Update the owner of a cellxgene collection.
    To run:
    ./scripts/cxg_admin.py --deployment prod update-collection-owner "$COLLECTION_ID $NEW_OWNER_ID
    """

    with db_session_manager() as session:
        key = (collection_id, CollectionVisibility.PUBLIC.name)
        collection = Collection.get(session, key)
        collection_name = collection.to_dict()["name"]

        if collection is not None:
            click.confirm(
                f"Are you sure you want to update the owner of the collection:{collection_name}?",
                abort=True,
            )
            collection.update(owner=new_owner)
            collection = Collection.get(session, key)
            if collection and (collection.owner == new_owner):
                click.echo(f"Updated owner of collection:{collection_id}. Owner is now {collection.owner}")
                exit(0)
            else:
                click.echo(f"Failed to update owner for collection_id:{collection_id}")
                exit(0)


def transfer_collections(ctx, curr_owner, new_owner):
    """Transfer all collections owned by the curr_owner to the new_owner.
    Retrieve user ids from auth0 before running or ping an engineer on the team to check the id of the owner in the
    database
    To run:
    ./scripts/cxg_admin.py --deployment prod transfer-collections $CURR_OWNER_ID $NEW_OWNER_ID
    """

    with db_session_manager() as session:
        collections = session.query(DbCollection).filter(DbCollection.owner == curr_owner).all()
        new_owner_collections_count = len(session.query(DbCollection).filter(DbCollection.owner == new_owner).all())
        if len(collections):
            click.confirm(
                f"Are you sure you want to update the owner of {len(collections)} collection"
                f"{'s' if len(collections) > 1 else ''} from {curr_owner} to "
                f"{new_owner}?",
                abort=True,
            )
            updated = (
                session.query(DbCollection)
                .filter(DbCollection.owner == curr_owner)
                .update({DbCollection.owner: new_owner})
            )
            session.commit()
            if updated > 0:
                collections = session.query(DbCollection).filter(DbCollection.owner == new_owner).all()
                click.echo(
                    f"{new_owner} previously owned {new_owner_collections_count}, they now own {len(collections)}"
                )
                click.echo(
                    f"Updated owner of collection for {updated} collections. {new_owner} is now the owner of "
                    f"{[[x.name, x.id] for x in collections]}"
                )
                exit(0)
            else:
                click.echo(
                    f"Failed to update owner for collections. {curr_owner} is still the owner of {len(collections)} "
                    f"collections"
                )
                exit(0)


def strip_all_collection_fields(ctx):
    """
    Strip all the `collection` string fields, so whitespace at the beginning and the end are removed.
    """
    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will strip whitespaces to all the string"
            "fields in the `collection` table",
            abort=True,
        )

        query = """
            update project
            set
            owner=TRIM(owner),
            name=TRIM(name),
            description=TRIM(description),
            contact_name=TRIM(contact_name),
            contact_email=TRIM(contact_email),
            data_submission_policy_version=TRIM(data_submission_policy_version)
            where
            owner != TRIM(owner) OR
            name != TRIM(name) OR
            description != TRIM(description) OR
            contact_name != TRIM(contact_name) OR
            contact_email != TRIM(contact_email) OR
            data_submission_policy_version != TRIM(data_submission_policy_version)
            """

        session.execute(query)
        session.commit()
