import os
from chalice import Chalice

from backend.corpora.lambdas.cloudfront_invalidator.cloudfront import (
    get_cloudfront_distribution,
    invalidate_distributions,
)


app = Chalice(app_name="cloudfront_invalidator")


@app.on_s3_event(
    name="corpora-frontend",
    bucket=os.getenv("S3_WEBSITE"),
    events=["s3:ObjectCreated:Post", "s3:ObjectCreated:Put"],
    prefix="index.html",
)
def handle_s3_event(event):
    print("EVENT-BUCKET:", event.bucket)
    distributions = get_cloudfront_distribution(event.bucket)
    print("DISTRIBUTIONS:", distributions)
    responses = invalidate_distributions(distributions)
    print("RESPONSES:", responses)
