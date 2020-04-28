"""
A chalice application is used to invalidate the cache of a cloudfront distribution for
an s3 bucket hosting a website. When a change is detected in the s3 bucket, the
AWS Cloudfront distribution for that bucket is invalidated. This allows the changers
made to a website to become immediately visible.
"""
import os
import sys
from chalice import Chalice

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from code.cloudfront import get_cloudfront_distribution, invalidate_distributions


app = Chalice(app_name="browser-frontend")


@app.on_s3_event(
    name="update-site",
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
