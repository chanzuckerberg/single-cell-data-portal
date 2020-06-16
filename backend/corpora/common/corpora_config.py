from .utils.secret_config import SecretConfig


class CorporaConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__('corpora/corpora', **kwargs)


class CorporaDbConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(component_name="corpora/corpora", secret_name="database", **kwargs)
