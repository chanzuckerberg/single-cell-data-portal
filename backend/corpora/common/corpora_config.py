from .utils.secret_config import SecretConfig


class CorporaDbConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(component_name="corpora/corpora", secret_name="database", **kwargs)
