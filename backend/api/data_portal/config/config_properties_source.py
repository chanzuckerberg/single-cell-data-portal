import json
import logging
import os
from abc import abstractmethod

from backend.api.data_portal.common.utils.aws import AwsSecret

logger = logging.getLogger(__name__)


class ConfigPropertiesSource:
    """
    Base class for obtaining configuration properties from a given source, as implemented by a subclass. Lazy loads
    and caches the configuration properties that are loaded.
    """

    def __init__(self):
        self._props = None

    @abstractmethod
    def load(self) -> dict:
        return {}

    def get_prop(self, name):
        if self._props is None:
            self._props = self.load()

        value = self._props.get(name)
        if value is not None:
            logger.debug(f"loaded config prop {name} from {self.__class__.__name__}")
        return value

    @staticmethod
    def raise_error(name):
        raise RuntimeError(name + " property not known")


class ComposedConfigPropertiesSource(ConfigPropertiesSource):
    def __init__(self, *config_properties_sources: ConfigPropertiesSource):
        super().__init__()
        self._config_properties_sources = config_properties_sources

    def load(self):
        pass

    def get_prop(self, name):
        for cps in self._config_properties_sources:
            if value := cps.get_prop(name):
                return value
        return None  # self.raise_error(name)


class DefaultConfigPropertiesSource(ConfigPropertiesSource):
    def __init__(self, **props):
        super().__init__()
        self._default_props = props

    def load(self):
        return self._default_props


class AwsSecretConfigPropertiesSource(ConfigPropertiesSource):
    def __init__(self, secret_path):
        self._secret_path = secret_path
        super().__init__()

    def load(self):
        secret = AwsSecret(self._secret_path)
        props = json.loads(secret.value)
        print(f"loaded config props from aws secret {self._secret_path}, keys={props.keys()}")
        return props


# TODO: Maybe we'll use the below someday. Or maybe not. I can't delete my own code.

# class EnvConfigPropertiesSource(ConfigPropertiesSource):
#     def load(self):
#         # convert all env var names to lower case for consistent lookup
#         env = dict([(k.lower(), v) for k, v in os.environ.items()])
#         return env
#
# class FileConfigPropertiesSource(ConfigPropertiesSource):
#     def __init__(self, file_path):
#         super().__init__()
#         self._file_path = file_path
#
#     def load(self):
#         with open(self._file_path, 'r') as f:
#             props = json.load(f)
#             assert props.isinstance(dict)
#             print(f'loaded config from {self._file_path}, keys={props.keys()}')
#             return props
