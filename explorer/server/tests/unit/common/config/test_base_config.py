import unittest

from server.common.config.app_config import AppConfig, flatten
from server.tests import FIXTURES_ROOT
from server.tests.unit.common.config import ConfigTests


class BaseConfigTest(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret", multi_dataset__dataroot=FIXTURES_ROOT)
        self.config.complete_config()

        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(config_file_name=self.config_file_name, **kwargs)
        config = AppConfig(file_name)
        return config

    def test_mapping_creation_returns_map_of_server_and_dataset_config(self):
        config = AppConfig()
        mapping = flatten(config.default_config)
        self.assertIsNotNone(mapping["server__app__verbose"])
        self.assertIsNotNone(mapping["default_dataset__presentation__max_categories"])

    def test_changes_from_default_returns_list_of_nondefault_config_values(self):
        config = self.get_config(verbose="true", lfc_cutoff=0.05)
        changes = config.changes_from_default()

        self.assertCountEqual(
            changes,
            [
                ("server__adaptor__cxg_adaptor__tiledb_ctx__py.init_buffer_bytes", 10485760, 536870912),
                ("server__app__verbose", True, False),
                ("server__app__flask_secret_key", "secret", None),
                ("server__data_locator__s3_region_name", "us-east-1", True),
                ("default_dataset__diffexp__lfc_cutoff", 0.05, 0.01),
                ("server__app__port", 5005, None),
                ("server__app__csp_directives", {}, None),
            ],
        )
