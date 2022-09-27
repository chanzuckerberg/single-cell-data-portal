import logging
from logging.config import dictConfig

from pythonjsonlogger import jsonlogger

from backend.corpora.common.logging_config import LOG_FORMAT, DATETIME_FORMAT


# The fields to log using the json logger.


def configure_logging(app_name):
    """https://docs.python.org/3/library/logging.config.html"""
    gunicorn_logger = logging.getLogger("gunicorn.error")
    dictConfig(
        {
            "version": 1,  # The version of dictConfig to use. This must be 1.
            "formatters": {
                "default": {"format": LOG_FORMAT, "()": jsonlogger.JsonFormatter, "datefmt": DATETIME_FORMAT}
            },
            "handlers": {
                "wsgi": {
                    "class": "logging.StreamHandler",
                    "stream": "ext://flask.logging.wsgi_errors_stream",
                    "formatter": "default",
                    "level": gunicorn_logger.level,
                }
            },
            "loggers": {
                app_name: {"level": gunicorn_logger.level, "handlers": ["wsgi"], "propagate": 0},
            },
            "root": {"level": gunicorn_logger.level, "handlers": ["wsgi"]},
        }
    )
