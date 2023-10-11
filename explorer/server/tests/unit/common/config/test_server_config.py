import json
import os
import unittest
from unittest.mock import patch
from urllib.parse import quote

from server.common.config.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.common.utils.utils import find_available_port
from server.tests import FIXTURES_ROOT, PROJECT_ROOT
from server.tests.unit.common.config import ConfigTests


class TestServerConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret", multi_dataset__dataroot=FIXTURES_ROOT)
        self.server_config = self.config.server

        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}", config_file_name=self.config_file_name, **kwargs
        )
        config = AppConfig(file_name)
        return config

    def test_init_raises_error_if_default_config_is_invalid(self):
        with self.assertRaises(ConfigurationError):
            self.get_config(port="not_valid")

    def test_handle_app__throws_error_if_port_doesnt_exist(self):
        with self.assertRaises(ConfigurationError):
            self.get_config(port=99999999)

    @patch("server.common.config.config_model.discover_s3_region_name")
    def test_handle_data_locator_works_for_default_types(self, mock_discover_region_name):
        mock_discover_region_name.return_value = None
        # Default config
        self.assertEqual(self.config.server__data_locator__s3_region_name, None)
        # hard coded
        config = self.get_config()
        self.assertEqual(config.server__data_locator__s3_region_name, "us-east-1")
        # incorrectly formatted
        dataroots = {
            "d1": {"base_url": "set1", "dataroot": "/path/to/set1_datasets/"},
            "d2": {"base_url": "set2/subdir", "dataroot": "s3://shouldnt/work"},
        }
        file_name = self.custom_app_config(
            dataroots=dataroots, config_file_name=self.config_file_name, data_locator_region_name="true"
        )
        with self.assertRaises(ConfigurationError):
            config = AppConfig(file_name)

    @patch("server.common.config.config_model.discover_s3_region_name")
    def test_handle_data_locator_can_read_from_dataroot(self, mock_discover_region_name):
        mock_discover_region_name.return_value = "us-west-2"
        dataroots = {
            "d1": {"base_url": "set1", "dataroot": "/path/to/set1_datasets/"},
            "d2": {"base_url": "set2/subdir", "dataroot": "s3://hosted-cellxgene-dev"},
        }
        file_name = self.custom_app_config(
            dataroots=dataroots, config_file_name=self.config_file_name, data_locator_region_name="true"
        )
        config = AppConfig(file_name)
        self.assertEqual(config.server__data_locator__s3_region_name, "us-west-2")
        mock_discover_region_name.assert_called_once_with("s3://hosted-cellxgene-dev")

    def test_handle_app__sets_web_base_url(self):
        config = self.get_config(web_base_url="anything.com")
        self.assertEqual(config.server__app__web_base_url, "anything.com")

    def test_handle_data_source__errors_when_passed_zero(self):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}",
            config_file_name="two_data_roots.yml",
        )

        file_name = self.custom_app_config(config_file_name="zero_roots.yml")
        config = AppConfig(file_name)
        with self.subTest("zero roots"), self.assertRaises(ConfigurationError):
            config.handle_data_source()

    @unittest.skip("skip when running in github action")
    def test_get_api_base_url_works(self):
        # test the api_base_url feature, and that it can contain a path
        config = AppConfig()
        backend_port = find_available_port("localhost", 10000)
        config.update_server_config(
            app_port=backend_port,
            app__flask_secret_key="secret",
            app__api_base_url=f"http://localhost:{backend_port}/additional/path",
            multi_dataset__dataroot=f"{PROJECT_ROOT}/example-dataset",
        )
        server = self.create_app(config)
        server.testing = True
        session = server.test_client()
        response = session.get("/additional/path/d/pbmc3k.cxg/api/v0.2/config")

        self.assertEqual(response.status_code, 200)
        data_config = json.loads(response.data)
        self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")

        # test the health check at the correct url
        response = session.get("/additional/path/health")
        assert json.loads(response.data)["status"] == "pass"

    def test_get_web_base_url_works(self):
        tests = [
            (dict(web_base_url="www.thisisawebsite.com"), "www.thisisawebsite.com"),
            (dict(web_base_url="local", port=5001), "http://localhost:5001"),
            (dict(web_base_url="www.thisisawebsite.com/"), "www.thisisawebsite.com"),
            (dict(api_base_url="www.api_base.com/"), "www.api_base.com"),
        ]
        for params, expected in tests:
            with self.subTest(params):
                config = self.get_config(**params)
                web_base_url = config.server__app__web_base_url
                self.assertEqual(web_base_url, expected)

    def test_multi_dataset_raises_error_for_illegal_routes(self):
        dataroot = f"{PROJECT_ROOT}/example-dataset"
        # test for illegal url_dataroots
        for illegal in ("../b", "!$*", "\\n", "", "(bad)"):
            with self.subTest(illegal), self.assertRaises(ConfigurationError):
                self.config.update_config(
                    server__multi_dataset__dataroots={"d": {"base_url": illegal, "dataroot": dataroot}}
                )

    def test_multidataset_works_for_legal_routes(self):
        # test for legal url_dataroots
        for legal in ("d", "this.is-okay_", "a/b"):
            self.config.update_server_config(
                multi_dataset__dataroots={"d": {"base_url": legal, "dataroot": f"{PROJECT_ROOT}/example-dataset"}}
            )

    @patch("server.app.app.render_template")
    def test_mulitdatasets_work_e2e(self, mock_render_template):
        try:
            os.symlink(FIXTURES_ROOT, f"{FIXTURES_ROOT}/set2")
            os.symlink(FIXTURES_ROOT, f"{FIXTURES_ROOT}/set3")
        except FileExistsError:
            pass

        mock_render_template.return_value = "something"
        # test that multi dataroots work end to end
        self.config.update_config(
            server__multi_dataset__dataroots=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}/set2", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}/set3", base_url="set3"),
            ),
            default__dataset__app__about_legal_tos="tos_default.html",
        )

        # Change this default to test if the dataroot overrides below work.
        self.config.update_default_dataset_config(default__dataset__app__about_legal_tos="tos_default.html")

        server = self.create_app(self.config)
        server.testing = True
        session = server.test_client()

        def _get_v03_url(url):  # TODO inline and do not use an API call to generate
            response = session.get(f"{url}/api/v0.3/s3_uri")
            s3_uri = quote(quote(response.json, safe=""), safe="")
            return f"/s3_uri/{s3_uri}/api/v0.3"

        v03_url = _get_v03_url("/set1/1/2/pbmc3k.cxg")
        response = session.get(f"{v03_url}/config")
        data_config = json.loads(response.data)
        assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"

        v03_url = _get_v03_url("/set2/pbmc3k.cxg")
        response = session.get(f"{v03_url}/config")

        data_config = json.loads(response.data)
        self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")

        v03_url = _get_v03_url("/set3/pbmc3k.cxg")
        response = session.get(f"{v03_url}/config")
        data_config = json.loads(response.data)
        self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")

        response = session.get("/health")
        self.assertEqual(json.loads(response.data)["status"], "pass")

        # access a dataset (no slash)
        with self.subTest("access a dataset without a trailing a slash"):
            response = session.get("/set2/pbmc3k.cxg")
            self.assertEqual(response.status_code, 308)

        # access a dataset (with slash)
        with self.subTest("access a dataset with a slash"):
            response = session.get("/set2/pbmc3k.cxg/")
            self.assertEqual(response.status_code, 200)

        # cleanup
        os.unlink(f"{FIXTURES_ROOT}/set2")
        os.unlink(f"{FIXTURES_ROOT}/set3")

    @patch("server.dataset.cxg_dataset.CxgDataset.set_tiledb_context")
    def test_handle_adaptor(self, mock_tiledb_context):
        custom_config = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}",
            cxg_tile_cache_size=10,
            cxg_tiledb_py_init_buffer_size=10,
        )
        config = AppConfig(custom_config)
        config.handle_adaptor()
        mock_tiledb_context.assert_called_once_with(
            {
                "sm.tile_cache_size": 10,
                "py.init_buffer_bytes": 10,
                "vfs.s3.region": "us-east-1",
            }
        )
