import os
from time import sleep

import boto3
from botocore.errorfactory import ClientError


class AwsSecret:
    AWS_SECRETS_MGR_SETTLE_TIME_SEC = 2

    """
    Wrapper for AWS secrets.
    Usage:
        secret = AwsSecret(name="my/component/secret")
        secret.update(value='{"foo":"bar"}')
        secret.value
        # -> '{"foo":"bar"}'
        secret.delete()
    Update handles create vs update and undeletion if necessary.
    """

    debug_logging = False

    def __init__(self, name):
        self._debug("AwsSecret.__init__({})".format(name))
        self.name = name
        self.secrets_mgr = boto3.client(
            service_name="secretsmanager", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None
        )
        self.secret_metadata = None
        self._load()

    @property
    def value(self):
        if not self.exists_in_aws:
            raise RuntimeError("No such secret: {}".format(self.name))
        if self.is_deleted:
            raise RuntimeError("Secret {} is deleted".format(self.name))
        response = self.secrets_mgr.get_secret_value(SecretId=self.arn)
        return response["SecretString"]

    @property
    def exists_in_aws(self):
        return self.secret_metadata is not None

    @property
    def is_deleted(self):
        return self.exists_in_aws and "DeletedDate" in self.secret_metadata

    @property
    def arn(self):
        return self.secret_metadata["ARN"]

    def update(self, value):
        if not self.exists_in_aws:
            self._debug("AwsSecret.update({}) creating...".format(self.name))
            self.secrets_mgr.create_secret(Name=self.name, SecretString=value)
            self._load()
        else:
            if self.is_deleted:
                self._debug("AwsSecret.update({}) restoring...".format(self.name))
                self._restore()
            self._debug("AwsSecret.update({}) updating...".format(self.name))
            self.secrets_mgr.put_secret_value(SecretId=self.arn, SecretString=value)

    def delete(self):
        if not self.exists_in_aws:
            raise RuntimeError("No such secret: {}".format(self.name))
        if not self.is_deleted:
            self.secrets_mgr.delete_secret(SecretId=self.name)
            sleep(self.AWS_SECRETS_MGR_SETTLE_TIME_SEC)  # eventual consistency
            self._load()

    def _load(self):
        try:
            response = self.secrets_mgr.describe_secret(SecretId=self.name)
            if response:
                self._debug("AwsSecret.load({}) it exists".format(self.name))
                self.secret_metadata = response
        except ClientError as e:
            if e.response["Error"]["Code"] == "ResourceNotFoundException":
                #  Normal operation, secret does not exist yet.
                self._debug("AwsSecret.load({}) ResourceNotFoundException".format(self.name))
            else:
                raise e

    def _restore(self):
        self.secrets_mgr.restore_secret(SecretId=self.arn)
        self._load()

    def _debug(self, message):
        if self.debug_logging:
            print(message)


s3 = boto3.resource("s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)


def delete_many_from_s3(bucket_name: str, prefix: str) -> None:
    """
    This deletes everything with a specific prefix from the given bucket

    """
    if not prefix:  # protect from deleting everything
        raise ValueError
    bucket = s3.Bucket(bucket_name)
    bucket.objects.filter(Prefix=f"{prefix}").delete()
