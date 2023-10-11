import logging
import unittest
from os import listdir, path

from flask_compress import Compress
from flask_cors import CORS

from server.app.app import Server
from server.common.config.app_config import AppConfig
from server.tests import FIXTURES_ROOT


def app_config(extra_server_config=None, extra_dataset_config=None):
    extra_server_config = extra_server_config or {}
    extra_dataset_config = extra_dataset_config or {}
    config = AppConfig()
    config.update_server_config(
        app__flask_secret_key="secret",
        limits__diffexp_cellcount_max=None,
        limits__column_request_max=None,
        multi_dataset__dataroot="/test/",
    )
    config.update_default_dataset_config(
        embeddings__names=["umap", "tsne", "pca"], presentation__max_categories=100, diffexp__lfc_cutoff=0.01
    )
    config.update_server_config(**extra_server_config)
    config.update_default_dataset_config(**extra_dataset_config)
    config.complete_config()
    return config


class TestServer(Server):
    def __init__(self, app_config: AppConfig):
        super().__init__(app_config)
        self._extract_base_url_and_dataset_for_api_calls(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        app.config["COMPRESS_MIMETYPES"] = [
            "text/html",
            "text/css",
            "text/xml",
            "application/json",
            "application/javascript",
            "application/octet-stream",
        ]
        Compress(app)
        if app_config.server__app__debug:
            CORS(app, supports_credentials=True, origins=[r"^http://localhost:\d+"])

    @staticmethod
    def _extract_base_url_and_dataset_for_api_calls(app_config: AppConfig):
        """
        Convenience method for generating the url for the developer to use to access the front end, which is
        dependent on the name of the dataset being used (and the CXG_CLIENT_PORT, if set): .../<base_url>/<dataset>
        @param app_config: the AppConfig
        @return: None
        """
        dataroot_and_base_url_pairs: list = list(app_config.server__multi_dataset__dataroots.values())
        if len(dataroot_and_base_url_pairs) > 1:
            logging.warning("Found more than one dataroot -- will use first")
        first_pair: dict = dataroot_and_base_url_pairs[0]
        base_url: str = first_pair["base_url"]
        dataroot: str = first_pair["dataroot"]
        logging.info(f"Using base_url {base_url} and dataroot {dataroot}")
        try:
            files: list = listdir(dataroot)
            if len(files) > 1:
                logging.warning(f"Found more than one dataset in {dataroot}")
            for dataroot_file in files:
                if dataroot_file.endswith(".cxg"):
                    cxg_file = dataroot_file
                    logging.info(f"Explorer will load dataset: {path.join(dataroot, cxg_file)}")
                    with open(".test_base_url.txt", "a") as test_base_url_file:
                        test_base_url_file.write(f"{base_url}/{cxg_file}")
                    break

            with open(".test_server_port.txt", "w") as dataroot_file:
                dataroot_file.write(f"{app_config.server__app__port}")
        except FileNotFoundError:
            logging.warning(f"Unable to access {dataroot}. Make sure your dataroot exists locally.")


class BaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls, app_config=None):
        cls.TEST_DATASET_URL_BASE = "/d/pbmc3k.cxg"
        cls.TEST_URL_BASE = f"{cls.TEST_DATASET_URL_BASE}/api/v0.2/"
        cls.maxDiff = None
        cls.app = cls.create_app(app_config)

    @classmethod
    def create_app(cls, app_config=None):
        if not app_config:
            app_config = AppConfig()
            app_config.update_server_config(
                app__flask_secret_key="testing",
                app__debug=True,
                multi_dataset__dataroot=f"{FIXTURES_ROOT}",
                multi_dataset__index=True,
            )
        app_config.complete_config()

        app = TestServer(app_config).app

        app.testing = True
        app.debug = True

        return app
