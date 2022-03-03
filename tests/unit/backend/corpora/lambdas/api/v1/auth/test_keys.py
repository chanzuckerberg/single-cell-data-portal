import unittest

import backend.corpora.common.auth0_management_session
from backend.corpora.lambdas.api.v1.auth import keys


class TestAuth0ManagementSession(unittest.TestCase):

    @unittest.mock.patch("backend.corpora.lambdas.api.v1.auth.get")
    def test_session_refresh(self):
        session = backend.corpora.common.auth0_management_session.Auth0ManagementSession()
        session.headers.update({"Authorization": f"Bearer deliberate-wrong-token"})
