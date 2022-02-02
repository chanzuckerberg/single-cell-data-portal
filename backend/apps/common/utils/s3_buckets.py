import os

import boto3


class _Buckets:
    _portal_resource = None

    @property
    def portal_resource(self):
        if not self._portal_resource:
            self._portal_resource = boto3.resource(
                "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
            )
        return self._portal_resource

    _portal_client = None

    @property
    def portal_client(self):
        if not self._portal_client:
            self._portal_client = boto3.client(
                "s3",
                endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None,
                config=boto3.session.Config(signature_version="s3v4"),
            )
        return self._portal_client

    _explorer_bucket = None

    @property
    def explorer_bucket(self):
        if not self._explorer_bucket:
            self._explorer_bucket = self.portal_resource.Bucket(
                os.getenv("CELLXGENE_BUCKET", f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}")
            )
        return self._explorer_bucket


buckets = _Buckets()
