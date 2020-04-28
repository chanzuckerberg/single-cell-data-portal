from datetime import datetime
from typing import List

import boto3

cloudfront_client = boto3.client("cloudfront")


def get_cloudfront_distribution(bucket_name: str) -> List[str]:
    cf_distributions = cloudfront_client.list_distributions()["DistributionList"]["Items"]
    distributions = []
    for d in cf_distributions:
        for o in d["Origins"]["Items"]:
            if bucket_name in o["DomainName"]:
                distributions.append(d["Id"])
    return distributions


def invalidate_distributions(distributions: List[str]) -> List[dict]:
    responses = []
    for d in distributions:
        timestamp = datetime.utcnow().replace(second=0, microsecond=0).strftime("%Y%m%d%H%M%S")
        # Limit the number of invalidation request to 1 per minute per distribution.
        responses.append(
            cloudfront_client.create_invalidation(
                DistributionId=d,
                InvalidationBatch={"Paths": {"Quantity": 1, "Items": ["/*"]}, "CallerReference": timestamp},
            )
        )
    return responses
