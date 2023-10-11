import json
import os
import tempfile
import unittest
from urllib.parse import quote

from server.common.config.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.tests import FIXTURES_ROOT, PROJECT_ROOT
from server.tests.unit.common.config import ConfigTests


class TestDatasetConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.app_config = AppConfig()
        self.app_config.update_server_config(app__flask_secret_key="secret", multi_dataset__dataroot=FIXTURES_ROOT)
        self.dataset_config = self.app_config.default_dataset
        self.app_config.complete_config()
        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}", config_file_name=self.config_file_name, **kwargs
        )
        app_config = AppConfig(file_name)
        return app_config

    def test_init_datatset_config_sets_vars_from_default_config(self):
        app_config = AppConfig()
        self.assertEqual(app_config.default_dataset__presentation__max_categories, 1000)
        self.assertEqual(app_config.default_dataset__diffexp__lfc_cutoff, 0.01)

    def test_app_sets_script_vars(self):
        app_config = self.get_config(scripts=["path/to/script"])

        self.assertEqual(app_config.default_dataset__app__scripts, [{"src": "path/to/script"}])

        app_config = self.get_config(scripts=[{"src": "path/to/script", "more": "different/script/path"}])
        self.assertEqual(
            app_config.default_dataset__app__scripts, [{"src": "path/to/script", "more": "different/script/path"}]
        )

        app_config = self.get_config(scripts=["path/to/script", "different/script/path"])
        # TODO @madison -- is this the desired functionality?
        self.assertEqual(
            app_config.default_dataset__app__scripts, [{"src": "path/to/script"}, {"src": "different/script/path"}]
        )

        with self.assertRaises(ConfigurationError):
            app_config = self.get_config(scripts=[{"more": "different/script/path"}])

    def test_multi_dataset(self):
        app_config = AppConfig()
        try:
            os.symlink(FIXTURES_ROOT, f"{FIXTURES_ROOT}/set2")
            os.symlink(FIXTURES_ROOT, f"{FIXTURES_ROOT}/set3")
        except FileExistsError:
            pass

        # test that multi dataroots work end to end
        app_config.update_server_config(
            app__flask_secret_key="secret",
            multi_dataset__dataroots=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}/set2", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}/set3", base_url="set3"),
            ),
        )

        # Change this default to test if the dataroot overrides below work.
        app_config.update_default_dataset_config(app__about_legal_tos="tos_default.html")
        server = self.create_app(app_config)

        server.testing = True
        session = server.test_client()

        def _get_v03_url(url):  # TODO inline and do not use an API call to generate
            response = session.get(f"{url}/api/v0.3/s3_uri")
            s3_uri = quote(quote(response.json, safe=""), safe="")
            return f"/s3_uri/{s3_uri}/api/v0.3"

        v03_url = _get_v03_url("/set1/1/2/pbmc3k.cxg")
        with self.subTest("Test config for dataroot /set1/1/2/ returns the s1 config"):
            response1 = session.get(f"{v03_url}/config")
            data_config_set_1 = json.loads(response1.data)
            self.assertEqual(data_config_set_1["config"]["displayNames"]["dataset"], "pbmc3k")

        v03_url = _get_v03_url("/set2/pbmc3k.cxg")
        with self.subTest("Test config for dataroot /set2 returns the s2 config"):
            response2 = session.get(f"{v03_url}/config")
            data_config_set_2 = json.loads(response2.data)
            self.assertEqual(data_config_set_2["config"]["displayNames"]["dataset"], "pbmc3k")

        v03_url = _get_v03_url("/set3/pbmc3k.cxg")
        with self.subTest("Test config for dataroot /set3/ returns the default dataset config"):
            response3 = session.get(f"{v03_url}/config")
            data_config_set_3 = json.loads(response3.data)
            self.assertEqual(data_config_set_3["config"]["displayNames"]["dataset"], "pbmc3k")
            self.assertEqual(data_config_set_3["config"]["parameters"]["about_legal_tos"], "tos_default.html")

        response = session.get("/health")
        self.assertEqual(json.loads(response.data)["status"], "pass")
        # cleanup
        os.unlink(f"{FIXTURES_ROOT}/set2")
        os.unlink(f"{FIXTURES_ROOT}/set3")

    def test_configfile_with_specialization(self):
        # test that per_dataset_config config load the default config, then the specialized config

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    multi_dataset:
                        dataroots:
                            test:
                                base_url: test
                                dataroot: fake_dataroot

                """
                fconfig.write(config)
            app_config = AppConfig(configfile)

            # test config from specialization
            self.assertEqual(app_config.server__multi_dataset__dataroots["test"]["base_url"], "test")

    def test_X_approximate_distribution(self):
        with self.subTest("OK"):
            self.app_config.update_default_dataset_config(X_approximate_distribution="count")

        tests = ["auto", "bad"]
        for test in tests:
            with self.subTest(test), self.assertRaises(ConfigurationError):
                self.app_config.update_default_dataset_config(X_approximate_distribution=test)
