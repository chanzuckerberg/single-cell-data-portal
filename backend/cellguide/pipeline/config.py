import os

from backend.common.utils.secret_config import SecretConfig


class CellGuideConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="cellguide_config", **kwargs)

    def get_defaults_template(self):
        deployment_stage = os.getenv("DEPLOYMENT_STAGE", "test")
        defaults_template = {"bucket": f"cellguide-{deployment_stage}"}
        return defaults_template
