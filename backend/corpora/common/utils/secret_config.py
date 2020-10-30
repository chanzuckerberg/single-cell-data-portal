import json
import os

from .aws_secret import AwsSecret


class SecretConfig:
    """
    This Config class stores its configuration in an AwsSecret.

    Subclass it to create your config, as follows (Python 3.x):
      class MyComponentConfig(Config):
          def __init__(self, *args, **kwargs):
              super().__init__('my_component', **kwargs)

    Setup your secrets in AWS Secrets Manager (e.g. use Terraform), as a JSON hash, at path:
        corpora/my_component/<deployment>/secrets

    To access secrets, (assuming that 'my_secret' is a key of your JSON hash):
        MyComponentConfig().my_secret

    Implements singleton-like behavior using techniques from
    http://python-3-patterns-idioms-test.readthedocs.io/en/latest/Singleton.html. All instances of this class will share
    the same config data. If you subclass this class, the subclass gets its own data, but all instances of the subclass
    share that data.
    """

    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, "_config"):
            cls.reset()
        return super(SecretConfig, cls).__new__(cls)

    def __init__(self, component_name, deployment=None, source=None, secret_name="secrets"):
        """
        If source is specified, it must be the path to a JSON file
        """

        super(SecretConfig, self).__init__()

        self._component_name = component_name
        self._deployment = deployment or os.environ["DEPLOYMENT_STAGE"]
        self._secret_name = secret_name
        self._source = self._determine_source(source)

    @property
    def config(self):
        return self.__class__._config

    def __getattr__(self, name):
        if self.config_is_loaded():
            return (
                self.value_from_config(name)
                or self.value_from_env(name)
                or self.value_from_defaults(name)
                or self.raise_error(name)
            )
        else:
            return (
                self.value_from_env(name)
                or (self.load() and self.value_from_config(name))
                or self.value_from_defaults(name)
                or self.raise_error(name)
            )

    @classmethod
    def reset(cls):
        cls._config = None
        cls._defaults = {}
        cls.use_env = False

    def set(self, config):
        """
        Bypass the load mechanism and set secrets directly.
        Used in testing.
        """

        self.__class__._config = config
        self.__class__.use_env = False
        self.update_defaults()

    def get_defaults_template(self):
        return {}

    def update_defaults(self):
        for k, v in self.get_defaults_template().items():
            self.__class__._defaults[k] = v.format(**self.config)

    def load(self):
        if self._source == "aws":
            self.load_from_aws()
        else:
            self.load_from_file(self._source)
        self.update_defaults()
        return True  # so we can be used in 'and' statements

    def config_is_loaded(self):
        return self.config is not None

    def load_from_aws(self):
        secret_path = f"corpora/{self._component_name}/{self._deployment}/{self._secret_name}"
        secret = AwsSecret(secret_path)
        self.from_json(secret.value)

    def load_from_file(self, config_file_path):
        with open(config_file_path, "r") as config_fp:
            self.from_json(config_fp.read())

    def from_json(self, config_json):
        self.__class__._config = json.loads(config_json)

    def _determine_source(self, source):
        if source:
            pass
        elif "CONFIG_SOURCE" in os.environ:
            source = os.environ["CONFIG_SOURCE"]
        else:
            source = "aws"
        return source

    def value_from_config(self, name):
        if name in self.config:
            return self.config[name]
        else:
            return None

    def value_from_defaults(self, name):
        return self._defaults.get(name)

    def value_from_env(self, name):
        if self.__class__.use_env and name.upper() in os.environ:
            return os.environ[name.upper()]
        else:
            return None

    def raise_error(self, name):
        raise RuntimeError(name + " is not in configuration")
