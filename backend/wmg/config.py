import os

from backend.corpora.common.utils.secret_config import SecretConfig


class WmgConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="wmg_config", **kwargs)

    def get_defaults_template(self):
        deployment_stage = os.getenv("DEPLOYMENT_STAGE", "test")
        defaults_template = {"bucket": f"wmg-{deployment_stage}", "tiledb_mem_gb": "16"}
        return defaults_template
