import logging
import os
import re
import types
import gevent.monkey

from ddtrace import tracer
from ddtrace import patch_all

from ddtrace.filters import FilterRequestsOnUrl


def configure_datadog_tracing(deployment_stage):
    """
    Configure Datadog APM tracing.
    See https://ddtrace.readthedocs.io/en/stable/basic_usage.html#patch-all
    """

    if os.getenv('DD_AGENT_HOST') and os.getenv('DD_TRACE_AGENT_PORT'):
        logging.info(f"Configuring Datadog APM tracing (ddtrace)")
        log = logging.getLogger("ddtrace.tracer")
        log.setLevel(logging.DEBUG)

        # next line may be redundant with DD_GEVENT_PATCH_ALL env var in .happy/terraform/modules/service/main.tf
        gevent.monkey.patch_all()
        tracer.configure(
                hostname=os.environ['DD_AGENT_HOST'],
                port=os.environ['DD_TRACE_AGENT_PORT'],
                # Filter out health check endpoint (index page: '/')
                # TODO: This doesn't seem to be working, and we'd rather filter on the http agent, to identify the
                #  requests as coming from the load balancer, regardless of the URL (this URL is also the Swagger APIs
                #  landing page. See modules/service/main.tf DD_APM_FILTER_TAGS_REJECT env var
                settings={
                    'FILTERS': [
                        FilterRequestsOnUrl([r'http://.*/$']),
                    ],
                }
        )
        patch_all()

        # enable Datadog profiling for development
        if deployment_stage not in ['staging', 'prod']:
            # noinspection PyPackageRequirements,PyUnresolvedReferences
            import ddtrace.profiling.auto


class DatadogTraced(type):
    """
    When specified as the __metaclass__ of a class, wraps all methods with tracer.wrap(), thereby instrumenting each
    method for Datadog APM tracing.
    """

    def __new__(mcs, name, bases, attrs):

        for attr_name, attr_value in attrs.items():
            if (isinstance(attr_value, types.FunctionType)) \
                    and not re.match("__.*__", attr_name):
                print(f"Tracing function {attr_name}")
                attrs[attr_name] = tracer.wrap()(attr_value)

        return super().__new__(mcs, name, bases, attrs)
