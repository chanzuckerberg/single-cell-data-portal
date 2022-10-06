import boto3
import uuid
from typing import List

from backend.corpora.common.corpora_config import CorporaCloudfrontConfig

import logging

client = boto3.client("cloudfront")


# Since Cloudfront is only used in deployed environments (dev, staging, prod),
# only trigger an invalidation if the distribution_id is defined in secrets manager.
# Otherwise this will be a no-op
def create_invalidation(paths: List[str]):
    try:
        distribution = CorporaCloudfrontConfig().distribution_id
    except RuntimeError:  # Will be raised if the attribute is not found (i.e. in rdev)
        logging.debug("No Cloudfront distribution found in secrets, will not invalidate")
        return None
    return _create_invalidation(distribution, paths)


def _create_invalidation(distribution: str, paths: List[str]):
    invalidation_id = str(uuid.uuid4())
    logging.info(f"Requesting invalidation {invalidation_id} for distribution {distribution}")
    response = client.create_invalidation(
        DistributionId=distribution,
        InvalidationBatch={
            "Paths": {
                "Quantity": len(paths),
                "Items": paths,
            },
            "CallerReference": invalidation_id,
        },
    )
    logging.info(response)


def create_invalidation_for_index_paths():
    return create_invalidation(["/dp/v1/datasets/index", "/dp/v1/collections/index"])
