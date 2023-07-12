import logging
import os
import sys

from click import Context

from backend.curation.api.v1.curation.collections.common import validate_uuid_else_forbidden

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import CollectionId


def tombstone_collection(ctx: Context, collection_id: str) -> None:
    """
    Tombstones the collection specified by uuid.
    :param ctx: command context
    :param uuid: ID that identifies the collection to tombstone
    """
    try:
        validate_uuid_else_forbidden(collection_id)
    except Exception:
        logging.error(f"{collection_id} is not a valid uuid")
        exit(1)
    business_logic: BusinessLogic = ctx.obj["business_logic"]
    collection = business_logic.get_canonical_collection(CollectionId(collection_id))
    if not collection:
        logging.error(f"Collection {collection_id} does not exist")
        exit(1)
    elif collection.tombstoned:
        logging.error(f"Collection {collection_id} is already tombstoned")
        exit(1)
    business_logic.tombstone_collection(CollectionId(collection_id))
    print(f"Successfully tombstoned Collection {collection_id}")


def tombstone_dataset(ctx, uuid):
    """
    Remove a dataset from Cellxgene. This will delete its artifacts/genesets and mark the dataset as tombstoned so
     it no longer shows up in the data portal.
     You must first SSH into the target deployment using `make db/tunnel` before running.
      ./scripts/cxg_admin.py --deployment staging tombstone-dataset "57cf1b53-af10-49e5-9a86-4bc70d0c92b6"

    """
    pass
