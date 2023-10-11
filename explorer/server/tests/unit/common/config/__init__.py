import os
import random
import shutil

from server.common.config.app_config import AppConfig, flatten
from server.tests import FIXTURES_ROOT
from server.tests.unit import BaseTest


class ConfigTests(BaseTest):
    tmp_fixtures_directory = os.path.join(FIXTURES_ROOT, "tmp_dir")

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.tmp_fixtures_directory)

    @classmethod
    def setUpClass(cls) -> None:
        os.makedirs(cls.tmp_fixtures_directory, exist_ok=True)

    def custom_server_config(
        self,
        verbose="false",
        debug="false",
        host="localhost",
        port="null",
        open_browser="false",
        force_https="false",
        flask_secret_key="secret",
        generate_cache_control_headers="false",
        server_timing_headers="false",
        csp_directives="null",
        api_base_url="null",
        web_base_url="null",
        insecure_test_environment="false",
        client_id="null",
        client_secret="null",
        jwt_decode_options="null",
        session_cookie="true",
        cookie="null",
        dataroot="null",
        dataroots=None,
        index="false",
        data_locator_region_name="us-east-1",
        data_locator_api_base="null",
        gene_info_api_base="null",
        cxg_tile_cache_size=8589934592,
        cxg_tiledb_py_init_buffer_size=10485760,
        column_request_max=32,
        diffexp_cellcount_max="null",
        config_file_name="server_config.yaml",
    ):
        dataroots = {} if dataroot is None else dataroots
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        server_config_outline_path = os.path.join(FIXTURES_ROOT, "czi_hosted_server_config_outline.py")
        with open(server_config_outline_path, "r") as config_skeleton:
            config = config_skeleton.read()
            server_config = eval(config)
        with open(configfile, "w") as server_config_file:
            server_config_file.write(server_config)
        return configfile

    def custom_app_config(
        self,
        verbose="false",
        debug="false",
        host="localhost",
        port="null",
        open_browser="false",
        force_https="false",
        flask_secret_key="secret",
        generate_cache_control_headers="false",
        server_timing_headers="false",
        csp_directives="null",
        api_base_url="null",
        web_base_url="null",
        client_id="null",
        client_secret="null",
        jwt_decode_options="null",
        session_cookie="true",
        cookie="null",
        dataroot="null",
        dataroots=None,
        index="false",
        data_locator_region_name="us-east-1",
        data_locator_api_base="null",
        gene_info_api_base="null",
        cxg_tile_cache_size=8589934592,
        cxg_tiledb_py_init_buffer_size=10485760,
        column_request_max=32,
        diffexp_cellcount_max="null",
        scripts=None,
        inline_scripts=None,
        about_legal_tos="null",
        about_legal_privacy="null",
        max_categories=1000,
        custom_colors="true",
        enable_users_annotations="true",
        annotation_type="local_file_csv",
        hosted_file_directory="null",
        local_file_csv_directory="null",
        local_file_csv_file="null",
        embedding_names=None,
        enable_difexp="true",
        lfc_cutoff=0.01,
        top_n=10,
        environment=None,
        aws_secrets_manager_region=None,
        aws_secrets_manager_secrets=None,
        X_approximate_distribution="normal",
        config_file_name="app_config.yml",
    ):
        aws_secrets_manager_secrets = [] if aws_secrets_manager_secrets is None else aws_secrets_manager_secrets
        embedding_names = [] if embedding_names is None else embedding_names
        inline_scripts = [] if inline_scripts is None else inline_scripts
        scripts = [] if scripts is None else scripts
        dataroots = {} if dataroots is None else dataroots

        random_num = random.randrange(999999)
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        server_config = self.custom_server_config(
            verbose=verbose,
            debug=debug,
            host=host,
            port=port,
            open_browser=open_browser,
            force_https=force_https,
            flask_secret_key=flask_secret_key,
            generate_cache_control_headers=generate_cache_control_headers,
            server_timing_headers=server_timing_headers,
            csp_directives=csp_directives,
            api_base_url=api_base_url,
            web_base_url=web_base_url,
            client_id=client_id,
            client_secret=client_secret,
            jwt_decode_options=jwt_decode_options,
            session_cookie=session_cookie,
            cookie=cookie,
            dataroot=dataroot,
            dataroots=dataroots,
            index=index,
            data_locator_region_name=data_locator_region_name,
            data_locator_api_base=data_locator_api_base,
            gene_info_api_base=gene_info_api_base,
            cxg_tile_cache_size=cxg_tile_cache_size,
            cxg_tiledb_py_init_buffer_size=cxg_tiledb_py_init_buffer_size,
            column_request_max=column_request_max,
            diffexp_cellcount_max=diffexp_cellcount_max,
            config_file_name=f"temp_server_config_{random_num}.yml",
        )
        dataset_config = self.custom_dataset_config(
            scripts=scripts,
            inline_scripts=inline_scripts,
            about_legal_tos=about_legal_tos,
            about_legal_privacy=about_legal_privacy,
            max_categories=max_categories,
            custom_colors=custom_colors,
            enable_users_annotations=enable_users_annotations,
            annotation_type=annotation_type,
            hosted_file_directory=hosted_file_directory,
            local_file_csv_directory=local_file_csv_directory,
            local_file_csv_file=local_file_csv_file,
            embedding_names=embedding_names,
            enable_difexp=enable_difexp,
            lfc_cutoff=lfc_cutoff,
            top_n=top_n,
            X_approximate_distribution=X_approximate_distribution,
            config_file_name=f"temp_dataset_config_{random_num}.yml",
        )

        with open(configfile, "w") as app_config_file:
            with open(server_config) as fp:
                app_config_file.write(fp.read())
            with open(dataset_config) as fp:
                app_config_file.write(fp.read())

        return configfile

    def custom_dataset_config(
        self,
        scripts=None,
        inline_scripts=None,
        about_legal_tos="null",
        about_legal_privacy="null",
        max_categories=1000,
        custom_colors="true",
        enable_users_annotations="true",
        annotation_type="local_file_csv",
        hosted_file_directory="null",
        local_file_csv_directory="null",
        local_file_csv_file="null",
        embedding_names=None,
        enable_difexp="true",
        lfc_cutoff=0.01,
        top_n=10,
        X_approximate_distribution="normal",
        config_file_name="dataset_config.yml",
    ):
        script = [] if scripts is None else scripts
        inline_scripts = [] if inline_scripts is None else inline_scripts
        embedding_names = [] if embedding_names is None else embedding_names
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        dataset_config_outline_path = os.path.join(FIXTURES_ROOT, "czi_hosted_dataset_config_outline.py")
        with open(dataset_config_outline_path, "r") as config_skeleton:
            config = config_skeleton.read()
            dataset_config = eval(config)
        with open(configfile, "w") as dataset_config_file:
            dataset_config_file.write(dataset_config)

        return configfile

    @staticmethod
    def compare_configs(a: AppConfig, b: AppConfig):
        _a = flatten(a.config)
        _b = flatten(b.config)

        diff = []
        keys = {*_a.keys(), *_b.keys()}
        for key in keys:
            value_a, value_b = _a.get(key), _b.get(key)
            if value_a != value_b:
                diff.append((key, value_a, value_b))
        return diff
