import uuid

import boto3
import os


class ExistingAwsSecretTestFixture:
    """
    A test fixture to manage an AWS Secrets Manager secret.
    It is a context manager that will instantiate and tear down the secret around your code.
    Use it with:
        with ExistingAwsSecretTestFixture(secret_name="foo", secret_value="bar") as secret:
            print("My name is {secret.name}")
            # test something
    If you leave the secret_name blank, a unique name will be generated for you.
    """

    SECRET_ID_TEMPLATE = "corpora/cicd/dev/secrets/{}"
    EXISTING_SECRET_DEFAULT_VALUE = '{"top":"secret"}'

    def __init__(self, secret_name=None, secret_value=None):
        self.secrets_mgr = boto3.client("secretsmanager", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"))
        self.name = secret_name or self.SECRET_ID_TEMPLATE.format(uuid.uuid4())
        self._value = secret_value or self.EXISTING_SECRET_DEFAULT_VALUE
        self._secret = None

    def __enter__(self):
        self._secret = self.secrets_mgr.create_secret(Name=self.name, SecretString=self._value)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.delete()
        except:
            print("failed to delete secret")

    @property
    def value(self):
        return self._value

    @property
    def arn(self):
        return self._secret["ARN"]

    def delete(self):
        self.secrets_mgr.delete_secret(SecretId=self.name, ForceDeleteWithoutRecovery=True)
