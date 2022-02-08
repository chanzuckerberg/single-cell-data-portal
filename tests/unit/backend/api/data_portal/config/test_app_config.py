import unittest

from backend.api.data_portal.config.app_config import AuthConfig


class TestAppConfig(unittest.TestCase):
    def test_auth_config_is_singleton(self):
        assert AuthConfig() is AuthConfig()

    def test__can_update_auth_config_prop(self):
        ac = AuthConfig()
        ac.api_base_url = "new_api_base_url"
        assert AuthConfig().api_base_url == "new_api_base_url"

    def test__when_auth_api_base_url_updated_then_derived_props_are_updated(self):
        ac = AuthConfig()

        ac.api_base_url = "updated_api_base_url"

        assert "updated_api_base_url" in ac.api_authorize_url
        assert "updated_api_base_url" in ac.api_token_url
        assert "updated_api_base_url" in ac.api_userinfo_url
        assert "updated_api_base_url" in ac.internal_url
        assert "updated_api_base_url" in "".join(ac.issuer)

    # def test__when_config_prop_is_from_env_then_prop_name_is_lowercased(self):
    #     os.environ['SOME_KEY'] = 'some_value'
    #
    #     ec = EnvConfigPropertiesSource()
    #
    #     assert ec.get_prop('some_key') == 'some_value'
    #     assert ec.get_prop('SOME_KEY') is None
    #
