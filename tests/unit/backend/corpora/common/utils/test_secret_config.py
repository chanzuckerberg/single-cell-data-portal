import unittest
import uuid

try:
    from unittest.mock import PropertyMock, patch
except ImportError:
    from mock import patch

import boto3

from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup, fixture_file_path
from tests.unit.backend.corpora.fixtures.existing_aws_secret_test_fixture import ExistingAwsSecretTestFixture

from backend.corpora.common.utils.secret_config import SecretConfig
from backend.corpora.common.utils.aws_secret import AwsSecret


class BogoComponentConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super(BogoComponentConfig, self).__init__("bogo_component", **kwargs)


class TestSecretConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # AwsSecret.debug_logging = True
        # To reduce eventual consistency issues, get everyone using the same Secrets Manager session
        cls.secrets_mgr = boto3.client("secretsmanager")
        cls.patcher = patch("backend.corpora.common.utils.aws_secret.boto3.client")
        boto3_client = cls.patcher.start()
        boto3_client.return_value = cls.secrets_mgr

    @classmethod
    def tearDownClass(cls):
        AwsSecret.debug_logging = False
        cls.patcher.stop()

    def setUp(self):
        self.deployment_env = "bogo_env_{}".format(uuid.uuid4())
        self.secret_name = ExistingAwsSecretTestFixture.SECRET_ID_TEMPLATE.format(
            f"bogo_component/{self.deployment_env}")
        BogoComponentConfig.reset()

    def tearDown(self):
        pass

    def test_from_file(self):
        with EnvironmentSetup({"CONFIG_SOURCE": fixture_file_path("bogo_config.js")}):
            config = BogoComponentConfig(deployment=self.deployment_env)
            self.assertEqual("value_from_file", config.secret1)

    def test_from_aws(self):
        with ExistingAwsSecretTestFixture(
                secret_name=self.secret_name, secret_value='{"secret1":"secret1_from_cloud"}'
        ):
            with EnvironmentSetup({"CONFIG_SOURCE": None}):
                config = BogoComponentConfig(deployment=self.deployment_env, source="aws")
                self.assertEqual("secret1_from_cloud", config.secret1)

    def test_custom_secret_name(self):
        custom_secret_name = ExistingAwsSecretTestFixture.SECRET_ID_TEMPLATE.format(
            f"bogo_component/{self.deployment_env}/custom-secret-name")
        with ExistingAwsSecretTestFixture(secret_name=custom_secret_name, secret_value='{"secret1":"custom"}'):
            class BogoComponentCustomConfig(SecretConfig):
                def __init__(self, *args, **kwargs):
                    super(BogoComponentCustomConfig, self).__init__(
                        "bogo_component", secret_name="custom-secret-name", **kwargs
                    )

            config = BogoComponentCustomConfig(deployment=self.deployment_env, source="aws")
            self.assertEqual("custom", config.secret1)

    @patch("backend.corpora.common.utils.aws_secret.AwsSecret")
    def test_singletonness(self, mock_AwsSecret):
        value_mock = PropertyMock(return_value='{"secret2": "foo"}')
        type(mock_AwsSecret()).value = value_mock

        config1 = BogoComponentConfig(deployment=self.deployment_env, source="aws")
        self.assertEqual("foo", config1.secret2)

        config2 = BogoComponentConfig(deployment=self.deployment_env, source="aws")
        self.assertEqual("foo", config2.secret2)

        value_mock.assert_called_once()

    # TRUTH TABLE
    # ITEM IS IN CONFIG | ITEM IS IN ENV | use_env IS SET | RESULT
    #        no         |       no       |       no       | exception
    #        no         |       no       |       yes      | exception
    #        no         |       yes      |       no       | exception
    #        no         |       yes      |       yes      | return env value
    #        yes        |       no       |       no       | return config value
    #        yes        |       no       |       yes      | return config value
    #        yes        |       yes      |       no       | return config value
    #        yes        |       yes      |       no       | return config value
    #        yes        |       yes      |       yes      | return env value

    def test_when_item_is_not_in_config_not_in_env_we_raise(self):
        with EnvironmentSetup({"CONFIG_SOURCE": None}):
            with self.assertRaises(RuntimeError):
                config = BogoComponentConfig(deployment=self.deployment_env)
                print(config.secret_that_we_never_put_into_config)

    def test_when_item_is_not_in_config_but_is_in_env_and_use_env_is_not_set_we_raise(self):
        with ExistingAwsSecretTestFixture(secret_name=self.secret_name, secret_value="{}"):
            with EnvironmentSetup({"CONFIG_SOURCE": None, "SECRET1": "secret1_from_env"}):
                with self.assertRaises(RuntimeError):
                    config = BogoComponentConfig(deployment=self.deployment_env)
                    print(config.secret1)

    def test_when_item_is_not_in_config_but_is_in_env_and_use_env_is_set_we_use_env(self):
        with ExistingAwsSecretTestFixture(secret_name=self.secret_name, secret_value="{}"):
            BogoComponentConfig.use_env = True
            with EnvironmentSetup({"CONFIG_SOURCE": None, "SECRET1": "secret1_from_env"}):
                config = BogoComponentConfig(deployment=self.deployment_env)
                self.assertEqual("secret1_from_env", config.secret1)

    def test_when_item_is_in_config_but_not_in_env_and_use_env_is_not_set_we_use_config(self):
        with ExistingAwsSecretTestFixture(
                secret_name=self.secret_name, secret_value='{"secret1":"secret1_from_cloud"}'
        ):
            with EnvironmentSetup({"CONFIG_SOURCE": None}):
                config = BogoComponentConfig(deployment=self.deployment_env)
                self.assertEqual("secret1_from_cloud", config.secret1)

    def test_when_item_is_in_config_but_not_in_env_and_use_env_is_set_we_use_config(self):
        with ExistingAwsSecretTestFixture(
                secret_name=self.secret_name, secret_value='{"secret1":"secret1_from_cloud"}'
        ):
            BogoComponentConfig.use_env = True
            with EnvironmentSetup({"CONFIG_SOURCE": None, "SECRET1": "secret1_from_env"}):
                config = BogoComponentConfig(deployment=self.deployment_env)
                self.assertEqual("secret1_from_env", config.secret1)

    def test_when_item_is_in_config_and_is_in_env_and_use_env_is_set_we_use_env(self):
        with ExistingAwsSecretTestFixture(
                secret_name=self.secret_name, secret_value='{"secret1":"secret1_from_cloud"}'
        ):
            BogoComponentConfig.use_env = True
            with EnvironmentSetup({"CONFIG_SOURCE": None, "SECRET1": "secret1_from_env"}):
                config = BogoComponentConfig(deployment=self.deployment_env)
                self.assertEqual("secret1_from_env", config.secret1)
