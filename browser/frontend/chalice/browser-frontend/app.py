"""
A chalice application is used to invalidate the cache of a cloudfront distribution for
an s3 bucket hosting a website. When a change is detected in the s3 bucket, the
AWS Cloudfront distribution for that bucket is invalidated. This allows the changers
made to a website to become immediately visible.
"""
from typing import List

import boto3
from chalice import Chalice

app = Chalice(app_name="browser-frontend")
_cloudfront = None


@app.on_s3_event(
    name="update-site",
    bucket="dcp-static-site-dev-699936264352",
    events=["s3:ObjectCreated:Post", "s3:ObjectCreated:Put"],
    prefix="index.html",
)
def handle_s3_event(event):
    print("EVENT-BUCKET:", event.bucket)
    cf = get_cloudfront_distribution(event.bucket)
    print("DISTRIBUTIONS:", cf)


def get_cloudfront_distribution(bucket_name: str) -> List[str]:
    cf_distributions = _cloudfront.list_distributions()["DistributionList"]["Items"]
    cf = []
    for d in cf_distributions:
        for o in d["Origins"]["Items"]:
            if bucket_name in o["DomainName"]:
                cf.append(d["Id"])
    return cf


def get_cloudfront_client():
    global _cloudfront
    if not _cloudfront:
        _cloudfront = boto3.client("cloudfront")
    return _cloudfront
