import os
import unittest
from time import sleep

import boto3
from mock import patch

from backend.corpora.common.utils.aws_secret import AwsSecret
from tests.unit.backend.corpora.fixtures.existing_aws_secret_test_fixture import ExistingAwsSecretTestFixture


class TestAwsSecret(unittest.TestCase):
    UNKNOWN_SECRET = ExistingAwsSecretTestFixture.SECRET_ID_TEMPLATE.format(
        "secret_that_does_not_exist"
    )  # Don't ever create this

    @classmethod
    def setUpClass(cls):
        # To reduce eventual consistency issues, get everyone using the same Secrets Manager session
        cls.secrets_mgr = boto3.client("secretsmanager", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"))
        cls.patcher = patch("backend.corpora.common.utils.aws_secret.boto3.client")
        boto3_client = cls.patcher.start()
        boto3_client.return_value = cls.secrets_mgr

    @classmethod
    def tearDownClass(cls):
        cls.patcher.stop()

    def test_init_of_unknown_secret_does_not_set_secret_metadata(self):
        secret = AwsSecret(name=self.UNKNOWN_SECRET)
        self.assertEqual(secret.secret_metadata, None)

    def test_init_of_existing_secret_retrieves_secret_metadata(self):
        with ExistingAwsSecretTestFixture() as existing_secret:
            secret = AwsSecret(name=existing_secret.name)
            self.assertIsNotNone(secret.secret_metadata)

    # Value Test Cases

    def test_value_of_unknown_secret_raises_exception(self):
        with ExistingAwsSecretTestFixture():
            secret = AwsSecret(name=self.UNKNOWN_SECRET)
            with self.assertRaisesRegex(RuntimeError, "No such"):
                secret.value  # noqa

    def test_value_of_existing_deleted_secret_raises_exception(self):
        with ExistingAwsSecretTestFixture() as existing_secret:
            secret = AwsSecret(name=existing_secret.name)
            secret.delete()
            with self.assertRaisesRegex(RuntimeError, "deleted"):
                x = secret.value  # noqa

    def test_value_of_existing_secret_returns_value(self):
        with ExistingAwsSecretTestFixture() as existing_secret:
            secret = AwsSecret(name=existing_secret.name)
            self.assertEqual(secret.value, existing_secret.value)

    # Update Test Cases

    def test_delete_of_unknown_secret_raises_exception(self):
        secret = AwsSecret(name=self.UNKNOWN_SECRET)
        with self.assertRaises(RuntimeError):
            secret.delete()

    def test_update_of_existing_secret_updates_secret(self):
        with ExistingAwsSecretTestFixture() as existing_secret:
            secret = AwsSecret(name=existing_secret.name)
            secret.update(value='{"foo":"bar"}')
            sleep(AwsSecret.AWS_SECRETS_MGR_SETTLE_TIME_SEC)
            self.assertEqual(
                self.secrets_mgr.get_secret_value(SecretId=existing_secret.arn)["SecretString"], '{"foo":"bar"}'
            )

    # Delete Test Cases

    def test_delete_of_existing_secret_deletes_secret(self):
        with ExistingAwsSecretTestFixture() as existing_secret:
            secret = AwsSecret(name=existing_secret.name)
            secret.delete()
            secret = self.secrets_mgr.describe_secret(SecretId=existing_secret.arn)
            self.assertIn("DeletedDate", secret)
