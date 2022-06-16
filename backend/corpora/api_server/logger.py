import logging
from logging.config import dictConfig

from pythonjsonlogger import jsonlogger


# The fields to log using the json logger.
LOGGED_FIELDS = ["levelname", "asctime", "name", "message"]
LOG_FORMAT = " ".join([f"%({field})" for field in LOGGED_FIELDS])


def configure_logging():
    """https://docs.python.org/3/library/logging.config.html"""
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
            "root": {"level": logging.INFO, "handlers": ["wsgi"]},
        }
    )
