import logging
from logging.config import dictConfig

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOGGED_FIELDS, format_log_string


def configure_logging(app_name):
    """https://docs.python.org/3/library/logging.config.html"""
    gunicorn_logger: logging.Logger = logging.getLogger("gunicorn.error")
    # The fields to log using the json logger.
    log_format = format_log_string(LOGGED_FIELDS + ["request_id"])
    dictConfig(
        {
            "version": 1,  # The version of dictConfig to use. This must be 1.
            "filters": {
                "request_id": {
                    "()": "backend.api_server.request_id.RequestIdFilter",
                },
            },
            "formatters": {
                "default": {"format": log_format, "()": jsonlogger.JsonFormatter, "datefmt": DATETIME_FORMAT}
            },
            "handlers": {
                "wsgi": {
                    "class": "logging.StreamHandler",
                    "stream": "ext://flask.logging.wsgi_errors_stream",
                    "formatter": "default",
                    "level": gunicorn_logger.level,
                    "filters": ["request_id"],
                }
            },
            "loggers": {
                app_name: {"level": gunicorn_logger.level, "handlers": ["wsgi"], "propagate": 0},
            },
            "root": {"level": gunicorn_logger.level, "handlers": ["wsgi"]},
        }
    )
