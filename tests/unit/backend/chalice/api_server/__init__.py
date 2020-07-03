import os

from tests.unit.backend.chalice import ChaliceTestHarness

os.environ["APP_NAME"] = "corpora-api"


class BaseAPITest:
    """
    Provide access to the a Chalice app hosting the Corpora API. All test for the Corpora API should inherit this class.
    """

    @classmethod
    def setUpClass(cls):
        if not cls.packaged:
            corpora_api_dir = os.path.join(os.environ["CORPORA_HOME"], "backend", "chalice", "api_server")
            cls.packaged = True
        cls.app = ChaliceTestHarness(corpora_api_dir)
