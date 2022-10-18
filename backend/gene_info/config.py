from backend.common.utils.secret_config import SecretConfig


class GeneInfoConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="gene_info_config", **kwargs)

    def get_defaults_template(self):
        defaults_template = {"ncbi_api_key": ""}
        return defaults_template
