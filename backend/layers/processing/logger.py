import logging

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT


def configure_logging():
    formatter = jsonlogger.JsonFormatter(LOG_FORMAT, DATETIME_FORMAT)
    logging.basicConfig(formatter=formatter)
