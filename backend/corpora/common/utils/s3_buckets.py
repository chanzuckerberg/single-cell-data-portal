import os

import boto3

from backend.corpora.common.utils.singleton import Singleton


class Buckets(metaclass=Singleton):
    _resource = None

    @property
    def resource(self):
        if not self._resource:
            self._resource = boto3.resource(
                "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
            )
        return self._resource

    _client = None

    @property
    def client(self):
        if not self._client:
            self._client = boto3.client(
                "s3",
                endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None,
                config=boto3.session.Config(signature_version="s3v4"),
            )
        return self._client

    _cxg_bucket = None

    @property
    def cxg_bucket(self):
        if not self._cxg_bucket:
            self._cxg_bucket = self.resource.Bucket(
                os.getenv("CELLXGENE_BUCKET", f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}")
            )
        return self._cxg_bucket


buckets = Buckets()
