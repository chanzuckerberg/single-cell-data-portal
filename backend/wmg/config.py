from backend.corpora.api_server.app import DEPLOYMENT_STAGE
from backend.corpora.common.utils.secret_config import SecretConfig


class WmgConfig(SecretConfig):

    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="wmg_config", **kwargs)

    def get_defaults_template(self):
        defaults_template = {
            "bucket": f"wmg-{DEPLOYMENT_STAGE}"
        }
        return defaults_template


