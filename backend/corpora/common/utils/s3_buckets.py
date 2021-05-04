import os

import boto3

s3_client = boto3.client(
    "s3",
    endpoint_url=os.getenv("BOTO_ENDPOINT_URL"),
    config=boto3.session.Config(signature_version="s3v4"),
)
_s3 = boto3.resource(
    "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
)
cxg_bucket = _s3.Bucket(os.getenv("CELLXGENE_BUCKET", f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}"))
