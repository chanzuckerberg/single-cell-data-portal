from backend.corpora.common.utils.secret_config import SecretConfig


class GeneInfoConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="gene_info_config", **kwargs)

    def get_defaults_template(self):
        defaults_template = {"ncbi_api_key": ""}
        return defaults_template

    # TODO: promote this impl to parent class, if new behavior works universally
    def __getattr__(self, name):
        # Environment variables intentionally override config file.
        if not self.config_is_loaded():
            self.load()
        if (value := self.value_from_env(name)) is not None:
            return value
        if (value := self.value_from_config(name)) is not None:
            return value
        if (value := self.value_from_defaults(name)) is not None:
            return value
        self.raise_error(name)
