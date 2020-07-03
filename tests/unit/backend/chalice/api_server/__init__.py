import os

from tests.unit.backend.chalice import ChaliceTestHarness, run

os.environ["APP_NAME"] = "corpora-api"



class BaseAPITest:
    """
    Provide access to the a Chalice app hosting the Corpora API. All test for the Corpora API should inherit this class.
    """

    packaged = False  # flag to only package once.

    @classmethod
    def setUpClass(cls):
        if not cls.packaged:
            corpora_api_dir = os.path.join(os.environ["CORPORA_HOME"], "backend", "chalice", "api_server")
            #  Packaging the chalice corpora api server. This is be slow but only needs to run once.
            # run(["make", "package", "-C", corpora_api_dir])
            cls.packaged = True
        cls.app = ChaliceTestHarness(corpora_api_dir)
