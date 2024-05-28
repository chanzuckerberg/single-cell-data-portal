import os

from backend.common.utils.secret_config import SecretConfig


class CensusCubeConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="wmg_config", **kwargs)

    def get_defaults_template(self):
        deployment_stage = os.getenv("DEPLOYMENT_STAGE", "test")
        defaults_template = {"bucket": f"wmg-{deployment_stage}", "data_path_prefix": "", "tiledb_config_overrides": {}}
        return defaults_template
