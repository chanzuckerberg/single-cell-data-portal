import os

from .utils.secret_config import SecretConfig


class CorporaConfig(SecretConfig):
    environ_source = "CORPORA_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="config", **kwargs)

    def get_defaults_template(self):
        template = {"upload_file_formats": ["h5ad"], "upload_max_file_size_gb": 30}
        upload_snf_arn = os.getenv("UPLOAD_SFN_ARN")
        if upload_snf_arn:
            template["upload_sfn_arn"] = upload_snf_arn
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
        # instead of {"database_uri": "postgresql://blah"} so we can add a suffix here
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
        super().__init__(
            component_name="backend",
            secret_name="auth0-secret",
            **kwargs,
        )

    def get_defaults_template(self):
        template = {
            "api_authorize_url": "{api_base_url}/authorize",
            "api_token_url": "{api_base_url}/oauth/token",
            "api_userinfo_url": "{api_base_url}/userinfo",
            "api_auth0_v2_url": "{api_signin_url}api/v2/",
            "internal_url": "{api_base_url}",
            "issuer": [],
        }
        template["issuer"].append(self.api_base_url + "/" if not self.api_base_url.endswith("/") else self.api_base_url)
        if self.config.get("api_signin_url"):
            # Adding the API sign in URL to the list of allow token issues. This allow the API to accept Auth token
            # generated for testing. Used in dev and staging.
            template["issuer"].append(self.api_signin_url)

        return template
