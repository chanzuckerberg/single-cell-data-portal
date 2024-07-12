import os

import pytest

skip_creation_on_prod = pytest.mark.skipif(
    os.environ["DEPLOYMENT_STAGE"] == "prod", reason="Do not make test collections public in prod"
)
skip_no_explorer_in_rdev = pytest.mark.skipif(
    os.environ["DEPLOYMENT_STAGE"] == "rdev", reason="Explorer is not available in rdev"
)
