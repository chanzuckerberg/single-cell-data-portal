import os

from backend.api.data_portal.config.config_properties_source import (
    ComposedConfigPropertiesSource,
    AwsSecretConfigPropertiesSource,
    DefaultConfigPropertiesSource,
)


# TODO: All app configuration is performed herein using `*Config` singleton classes. Would it be simpler to just use
# globals variables instead of singleton classes? We already have an `app` global var (in app.py), and globals are
# a reasonably Python idiom to use when singleton-like behavior is required.


class ProcessingConfigPropertiesSource(ComposedConfigPropertiesSource):
    def __init__(self):
        deployment = os.environ["DEPLOYMENT_STAGE"]
        super().__init__(
            AwsSecretConfigPropertiesSource(f"corpora/backend/{deployment}/config"),
            DefaultConfigPropertiesSource(slack_webhook="", upload_max_file_size_gb=30, upload_file_formats=["h5ad"]),
        )


class ProcessingConfig:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(ProcessingConfig, cls).__new__(cls)
            cls._instance.__singleton_init()
        return cls._instance

    def __singleton_init(self):
        config_properties_source = ProcessingConfigPropertiesSource()

        self.__dict__["upload_max_file_size_gb"] = config_properties_source.get_prop("upload_max_file_size_gb")
        self.__dict__["upload_file_formats"] = config_properties_source.get_prop("upload_file_formats")
        self.__dict__["upload_sfn_arn"] = config_properties_source.get_prop("upload_sfn_arn")
        self.__dict__["slack_webhook"] = config_properties_source.get_prop("slack_webhook")
        # TODO: This prop is not exactly "processing"-related, as it is used by endpoints that serve the frontend
        self.__dict__["crossref_api_uri"] = config_properties_source.get_prop("crossref_api_uri")


class DbConfigPropertiesSource(AwsSecretConfigPropertiesSource):
    def __init__(self):
        deployment = os.environ["DEPLOYMENT_STAGE"]
        # HACK: Why have a _local_ config prop in AWS? This seems wrong. Eliminate?
        secret_name = f"database{'_local' if 'CORPORA_LOCAL_DEV' in os.environ else ''}"

        super().__init__(f"corpora/backend/{deployment}/{secret_name}")


class DbConfig:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(DbConfig, cls).__new__(cls)
            cls._instance.__singleton_init()
        return cls._instance

    def __singleton_init(self):
        config_properties_source = DbConfigPropertiesSource()
        self.__dict__["database_uri"] = config_properties_source.get_prop(
            "database_uri"
        ) or config_properties_source.get_prop(
            "remote_dev_uri"
        )  # HACK for rdev envs!


class AuthConfigPropertiesSource(AwsSecretConfigPropertiesSource):
    def __init__(self):
        deployment = os.environ["DEPLOYMENT_STAGE"]
        super().__init__(f"corpora/backend/{deployment}/auth0-secret")


class AuthConfig:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(AuthConfig, cls).__new__(cls)
            cls._instance.__singleton_init()
        return cls._instance

    def __singleton_init(self):
        config_properties_source = AuthConfigPropertiesSource()

        self._api_signin_url = config_properties_source.get_prop("api_signin_url")

        for source_prop in [
            "audience",
            "api_audience",
            "issuer",
            "cookie_name",
            "code_challenge_method",
            "redirect_to_frontend",
            "callback_base_url",
            "client_id",
            "client_secret",
            "api_base_url",
            "refresh_token_url",
            "access_token_url",
            "test_account_password",
            "test_account_username",
        ]:
            self.__dict__[source_prop] = config_properties_source.get_prop(source_prop)

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if key == "api_base_url":
            self.__update_api_base_url(value)

    def __update_api_base_url(self, api_base_url):
        """
        Update config props that are derived from the "api_base_url" property. This is only necessary to support
        updates while running tests, as some tests need to alter the api_base_url prop to perform a fake auth process
        """
        self.__dict__["api_authorize_url"] = f"{api_base_url}/authorize"
        self.__dict__["api_token_url"] = f"{api_base_url}/oauth/token"
        self.__dict__["api_userinfo_url"] = f"{api_base_url}/userinfo"
        self.__dict__["internal_url"] = api_base_url
        self.__dict__["issuer"] = self.build_issuers(api_base_url, self._api_signin_url)

    @staticmethod
    def build_issuers(api_base_url, api_signin_url):
        issuers = []
        if not api_base_url.endswith("/"):
            issuers.append(api_base_url + "/")
        else:
            issuers.append(api_base_url)

        if api_signin_url:
            # Adding the API sign in URL to the list of allow token issues. This allow the API to accept Auth token
            # generated for testing. Used in dev and staging.
            issuers.append(api_signin_url)

        return issuers


# class AppConfig:
#     """
#     Singleton app configuration
#     """
#
#     _instance = None
#
#     def __new__(cls):
#         if cls._instance is None:
#             cls._instance = app_config = super(AppConfig, cls).__new__(cls)
#             app_config._corpora_config = ProcessingConfig()
#             app_config._auth_config = AuthConfig()
#             app_config._db_config = DbConfig()
#         return cls._instance
#
#     @property
#     def corpora(self):
#         return self._corpora_config
#
#     @property
#     def auth(self):
#         return self._auth_config
#
#     @property
#     def db(self):
#         return self._db_config
