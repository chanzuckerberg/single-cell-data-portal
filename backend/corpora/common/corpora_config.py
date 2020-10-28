import os

from .utils.secret_config import SecretConfig


class CorporaConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("corpora/corpora", **kwargs)


class CorporaDbConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(
            component_name="backend",
            secret_name=f"database{'_local' if 'CORPORA_LOCAL_DEV' in os.environ else ''}",
            **kwargs,
        )


class CorporaAuthConfig(SecretConfig):
    """
    For a description of the secret key contents, see backend/config/auth0-secret-template.json.
    """

    def __init__(self, *args, **kwargs):
        deployment = os.environ["DEPLOYMENT_STAGE"]
        if deployment == "test":
            super().__init__(component_name="backend", deployment="dev", secret_name="auth0-secret", **kwargs)
            if not self.config_is_loaded():
                self.load()
            self.config["callback_base_url"] = "http://localhost:5000"
        else:
            super().__init__(
                component_name="backend",
                secret_name="auth0-secret",
                **kwargs,
            )
            if not self.config_is_loaded():
                self.load()

        self.set_defaults(
            {
                "api_authorize_url": f"{self.config['api_base_url']}/authorize",
                "api_token_url": f"{self.config['api_base_url']}/oauth/token",
            }
        )
