from datetime import datetime
from typing import List

import boto3

cloudfront_client = boto3.client("cloudfront")


def get_cloudfront_distribution(bucket_name: str) -> List[str]:
    cf_distributions = cloudfront_client.list_distributions()["DistributionList"]["Items"]
    distributions = []
    for distr in cf_distributions:
        for origin in distr["Origins"]["Items"]:
            if bucket_name in origin["DomainName"]:
                distributions.append(distr["Id"])
    return distributions


def invalidate_distributions(distributions: List[str]) -> List[dict]:
    """
    Using the AWS Cloudfront API to invalidate the cache of the deployed website. AWS requires all CallerReferences
    to be unique per invalidation request, otherwise an error is returned. This is to prevent identical requests see:
    https://docs.aws.amazon.com/cloudfront/latest/APIReference/API_InvalidationBatch.html. We use a time stamp to
    produce unique CallReferences and round the time stamp to the closest minute. This allows us to control how
    frequently an invalidation request can be made.
    :param distributions: The cloudfront distribution ID to invalidate
    :return: The responses from the invalidation requests.
    """

    responses = []
    for distr in distributions:
        timestamp = datetime.utcnow().replace(second=0, microsecond=0).strftime("%Y%m%d%H%M%S")
        # Produce a unique CallReference

        responses.append(
            cloudfront_client.create_invalidation(
                DistributionId=distr,
                InvalidationBatch={"Paths": {"Quantity": 1, "Items": ["/*"]}, "CallerReference": timestamp},
            )
        )
    return responses
