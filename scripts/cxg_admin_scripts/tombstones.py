import os
import sys

from click import Context

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import CollectionId


def tombstone_collection(ctx: Context, uuid: str) -> None:
    """
    Tombstones the collection specified by uuid.
    :param ctx: command context
    :param uuid: ID that identifies the collection to tombstone
    """
    business_logic: BusinessLogic = ctx.obj["business_logic"]
    business_logic.tombstone_collection(CollectionId(uuid))
    print(f"Successfully tombstoned Collection {uuid}")


def tombstone_dataset(ctx, uuid):
    """
    Remove a dataset from Cellxgene. This will delete its artifacts/genesets and mark the dataset as tombstoned so
     it no longer shows up in the data portal.
     You must first SSH into the target deployment using `make db/tunnel` before running.
      ./scripts/cxg_admin.py --deployment staging tombstone-dataset "57cf1b53-af10-49e5-9a86-4bc70d0c92b6"

    """
    pass
