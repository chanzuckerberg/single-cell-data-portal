import os
import tempfile
import unittest

import yaml

from server.common.config.app_config import AppConfig
from server.default_config import default_config
from server.tests import FIXTURES_ROOT
from server.tests.unit.common.config import ConfigTests


class AppConfigTest(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret", multi_dataset__dataroot=FIXTURES_ROOT)
        self.config.complete_config()
        self.maxDiff = None

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}", config_file_name=self.config_file_name, **kwargs
        )
        config = AppConfig(file_name)
        return config

    def test_get_default_config_correctly_reads_default_config_file(self):
        app_default_config = AppConfig().default_config

        expected_config = yaml.load(default_config, Loader=yaml.Loader)

        server_config = app_default_config["server"]
        dataset_config = app_default_config["default_dataset"]

        expected_server_config = expected_config["server"]
        expected_dataset_config = expected_config["default_dataset"]

        self.assertDictEqual(app_default_config, expected_config)
        self.assertDictEqual(server_config, expected_server_config)
        self.assertDictEqual(dataset_config, expected_dataset_config)

    def test_update_app_config(self):
        default_config = AppConfig()
        config = AppConfig()
        config.update_config(server__app__verbose=True, server__multi_dataset__dataroot="datadir")
        vars = self.compare_configs(config, default_config)
        self.assertCountEqual(
            vars,
            [
                ("server__app__verbose", True, False),
                ("server__multi_dataset__dataroots__d__base_url", "d", None),
                ("server__multi_dataset__dataroots__d__dataroot", "datadir", None),
            ],
        )

        config = AppConfig()
        config.update_config(default_dataset__app__scripts=(), default__dataset__app__inline_scripts=())
        vars = self.compare_configs(config, default_config)
        self.assertCountEqual(vars, [])

        config = AppConfig()
        config.update_config(default_dataset__app__scripts=("a", "b"), default_dataset__app__inline_scripts=["c", "d"])
        vars = self.compare_configs(config, default_config)
        self.assertCountEqual(
            vars,
            [
                ("default_dataset__app__scripts", [{"src": "a"}, {"src": "b"}], []),
                ("default_dataset__app__inline_scripts", ["c", "d"], []),
            ],
        )

    def test_configfile_no_dataset_section(self):
        # test a config file without a dataset section
        default_config = AppConfig()
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    app:
                        flask_secret_key: secret
                    multi_dataset:
                        dataroot: test_dataroot

                """
                fconfig.write(config)

            app_config = AppConfig(configfile)
            server_changes = self.compare_configs(app_config, default_config)
            self.assertCountEqual(
                server_changes,
                [
                    ("server__multi_dataset__dataroots__d__dataroot", "test_dataroot", None),
                    ("server__app__flask_secret_key", "secret", None),
                    ("server__multi_dataset__dataroots__d__base_url", "d", None),
                ],
            )

    def test_configfile_no_server_section(self):
        default_config = AppConfig()
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                default_dataset:
                  app:
                    about_legal_tos: expected_value
                """
                fconfig.write(config)

            app_config = AppConfig(configfile)
            changes = self.compare_configs(app_config, default_config)
            self.assertCountEqual(changes, [("default_dataset__app__about_legal_tos", "expected_value", None)])

    def test_csp_directives(self):
        default_config = AppConfig()
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                    server:
                      app:
                        csp_directives:
                          img-src:
                            - test_list
                          script-src: test_string
                """
                fconfig.write(config)

            app_config = AppConfig(configfile)
            changes = self.compare_configs(app_config, default_config)
            self.assertCountEqual(
                changes,
                [
                    ("server__app__csp_directives__img-src", ["test_list"], None),
                    ("server__app__csp_directives__script-src", ["test_string"], None),
                    ("server__app__csp_directives__connect-src", [], None),
                ],
            )
