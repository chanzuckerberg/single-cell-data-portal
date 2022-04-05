import boto3
import uuid

from backend.corpora.common.corpora_config import CorporaCloudfrontConfig

client = boto3.client('cloudfront')

# Since Cloudfront is only used in deployed environments (dev, staging, prod),
# only trigger an invalidation if the distribution_id is defined in secrets manager.
# Otherwise this will be a no-op
def create_invalidation(paths: list[str]):
    distribution = CorporaCloudfrontConfig().distribution_id
    if distribution:
        return _create_invalidation(distribution, paths)
    

def _create_invalidation(distribution: str, paths: list[str]):
    return client.create_invalidation(
        DistributionId=distribution,
        InvalidationBatch={
            'Paths': {
                'Quantity': len(paths),
                'Items': paths,
                },
            'CallerReference': str(uuid.uuid4())
            }
        )

def create_invalidation_for_index_paths():
    return create_invalidation(["/dp/v1/datasets/index", "/dp/v1/collections/index"])