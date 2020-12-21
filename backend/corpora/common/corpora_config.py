import os

from .utils.secret_config import SecretConfig


class CorporaConfig(SecretConfig):
    environ_source = "CORPORA_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="config", **kwargs)

    def get_defaults_template(self):
        template = {
            "upload_file_formats": ["h5ad"],
            "upload_max_file_size_gb": 30
        }
        if os.getenv("UPLOAD_SFN_ARN"):
            template["upload_sfn_arn"] = os.getenv("UPLOAD_SFN_ARN")
        return template


class CorporaDbConfig(SecretConfig):
    environ_source = "CORPORA_DB_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__(
            component_name="backend",
            secret_name=f"database{'_local' if 'CORPORA_LOCAL_DEV' in os.environ else ''}",
            **kwargs,
        )

    def get_defaults_template(self):
        # The db secret for remote dev envs is {"remote_dev_uri": "postgresql://blah"}
        # instead of {"database_url": "postgresql://blah"} so we can add a suffix here
        # based on the remote dev env name.
        remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX", "")
        if not remote_dev_prefix:
            return {}
        return {
            "database_uri": "{remote_dev_uri}" + remote_dev_prefix,
        }


class CorporaAuthConfig(SecretConfig):
    """
    For a description of the secret key contents, see backend/config/auth0-secret-template.json.
    """

    environ_source = "CORPORA_AUTH_CONFIG"

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

    def get_defaults_template(self):
        return {
            "api_authorize_url": "{api_base_url}/authorize",
            "api_token_url": "{api_base_url}/oauth/token",
            "internal_url": "{api_base_url}",
        }
