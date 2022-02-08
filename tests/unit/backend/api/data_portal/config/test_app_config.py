import unittest
from unittest.mock import patch

from backend.api.data_portal.config.app_config import AuthConfig, DbConfig


class TestAppConfig(unittest.TestCase):
    def test_auth_config_is_singleton(self):
        assert AuthConfig() is AuthConfig()

    def test__can_update_auth_config_prop(self):
        ac = AuthConfig()
        ac.api_base_url = "new_api_base_url"
        assert AuthConfig().api_base_url == "new_api_base_url"

    def test__when_auth_api_base_url_updated__then_derived_props_are_updated(self):
        ac = AuthConfig()

        ac.api_base_url = "updated_api_base_url"

        assert "updated_api_base_url" in ac.api_authorize_url
        assert "updated_api_base_url" in ac.api_token_url
        assert "updated_api_base_url" in ac.api_userinfo_url
        assert "updated_api_base_url" in ac.internal_url
        assert "updated_api_base_url" in "".join(ac.issuer)

    @patch("backend.api.data_portal.config.app_config.DbConfigPropertiesSource.load")
    @patch("backend.api.data_portal.config.app_config.os.getenv")
    def test__when_in_rdev_env__then_dbconfig_uses_rdev_db_uri(self, mock_os_get_env, mock_prop_source_load):
        mock_prop_source_load.return_value = dict(remote_dev_uri="some_rdev_db_uri_base/")
        mock_os_get_env.return_value = "rdev_env_name"

        DbConfig()._instance = None  # force re-init of this singleton obj

        assert DbConfig().database_uri == "some_rdev_db_uri_base/rdev_env_name"

    # def test__when_config_prop_is_from_env_then_prop_name_is_lowercased(self):
    #     os.environ['SOME_KEY'] = 'some_value'
    #
    #     ec = EnvConfigPropertiesSource()
    #
    #     assert ec.get_prop('some_key') == 'some_value'
    #     assert ec.get_prop('SOME_KEY') is None
    #
