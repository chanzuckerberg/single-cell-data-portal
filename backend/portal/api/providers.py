from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.batch_job_provider import BatchJobProvider
from backend.layers.thirdparty.cloudfront_provider import CloudfrontProvider
from backend.layers.thirdparty.crossref_provider import CrossrefProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProvider

_business_logic = None
_cloudfront_provider = None


def get_business_logic():
    global _business_logic
    if not _business_logic:
        _business_logic = BusinessLogic(
            DatabaseProvider(),
            BatchJobProvider(),
            CrossrefProvider(),
            StepFunctionProvider(),
            S3Provider(),
            UriProvider(),
        )
    return _business_logic


def get_cloudfront_provider():
    global _cloudfront_provider
    if not _cloudfront_provider:
        _cloudfront_provider = CloudfrontProvider()
    return _cloudfront_provider
