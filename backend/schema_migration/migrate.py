import os
import sys

from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider


def migrate(step_name: str, business_logic: BusinessLogic) -> bool:
    """
    Gets called by the step function at every different step, as defined by `step_name`
    """
    return True


if __name__ == "__main__":
    database_provider = DatabaseProvider()
    s3_provider = S3Provider()
    uri_provider = UriProvider()

    business_logic = BusinessLogic(
        database_provider,
        None,  # Not required
        None,  # Not required
        s3_provider,
        uri_provider,
    )
    step_name = os.environ["STEP_NAME"]
    rv = migrate(step_name, business_logic)
    sys.exit(rv)
