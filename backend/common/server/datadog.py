import os

from ddtrace import patch_all, tracer
from ddtrace.filters import FilterRequestsOnUrl

DEPLOYMENT_STAGE = os.environ["DEPLOYMENT_STAGE"]


def _should_configure_datadog_tracing():
    return (
        DEPLOYMENT_STAGE in ["dev", "staging", "prod"]
        and os.environ.get("DD_AGENT_HOST", None)
        and os.environ.get("DD_TRACE_AGENT_PORT", None)
    )


def initialize_datadog_tracing():
    if _should_configure_datadog_tracing():
        # Datadog APM tracing
        # See https://ddtrace.readthedocs.io/en/stable/basic_usage.html#patch-all

        tracer.configure(
            hostname=os.environ["DD_AGENT_HOST"],
            port=os.environ["DD_TRACE_AGENT_PORT"],
            # Filter out health check endpoint (index page: '/')
            settings={
                "FILTERS": [
                    FilterRequestsOnUrl([r"http://.*/$"]),
                ],
            },
        )
        patch_all()
