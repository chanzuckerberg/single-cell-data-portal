import unittest
from backend.corpora.lambdas.api.v1.curator import keys


class TestAuth0ManagementSession(unittest.TestCase):

    @unittest.mock.patch("backend.corpora.lambdas.api.v1.curator.get")
    def test_session_refresh(self):
        session = keys.Auth0ManagementSession()
        session.headers.update({"Authorization": f"Bearer deliberate-wrong-token"})
