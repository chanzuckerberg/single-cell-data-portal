import logging
from logging.config import dictConfig

from pythonjsonlogger import jsonlogger


# The fields to log using the json logger.
LOGGED_FIELDS = ["levelname", "asctime", "name", "message", "lineno", "module", "funcName", "pathname"]
LOG_FORMAT = " ".join([f"%({field})" for field in LOGGED_FIELDS])


def configure_logging(app_name):
    """https://docs.python.org/3/library/logging.config.html"""
    gunicorn_logger = logging.getLogger("gunicorn.error")
    dictConfig(
        {
            "version": 1,  # The version of dictConfig to use. This must be 1.
            "formatters": {
                "default": {"format": LOG_FORMAT, "()": jsonlogger.JsonFormatter, "datefmt": "%Y-%m-%dT%H:%M:%S.%03dZ"}
            },
            "handlers": {
                "wsgi": {
                    "class": "logging.StreamHandler",
                    "stream": "ext://flask.logging.wsgi_errors_stream",
                    "formatter": "default",
                    "level": "INFO",
                }
            },
            "loggers": {
                app_name: {"level": gunicorn_logger.level, "handlers": ["wsgi"]},
            },
            "root": {"level": gunicorn_logger.level, "handlers": ["wsgi"]},
        }
    )
