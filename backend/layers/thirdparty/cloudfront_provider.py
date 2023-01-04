import boto3
import uuid
from typing import List

from backend.common.corpora_config import CorporaCloudfrontConfig

import logging

from backend.layers.thirdparty.cdn_provider_interface import CDNProviderInterface

client = boto3.client("cloudfront")


class CloudfrontProvider(CDNProviderInterface):

    def create_invalidation(self, paths: List[str]):
        """
        Creates an invalidation for the specified paths.
        Since Cloudfront is only used in deployed environments (dev, staging, prod),
        only trigger an invalidation if the distribution_id is defined in secrets manager.
        Otherwise this will be a no-op
        """
        try:
            distribution = CorporaCloudfrontConfig().distribution_id
        except RuntimeError:  # Will be raised if the attribute is not found (i.e. in rdev)
            logging.debug("No Cloudfront distribution found in secrets, will not invalidate")
            return None
        return self._create_invalidation(distribution, paths)

    def _create_invalidation(self, distribution: str, paths: List[str]):
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

    def create_invalidation_for_index_paths(self):
        """
        Creates an invalidation for the index paths.
        By default, these are the only portal paths that require caching.
        """
        return self.create_invalidation(["/dp/v1/datasets/index", "/dp/v1/collections/index"])
