import logging
import sys

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT


def configure_logging(level=logging.INFO):
    log_stdout_handler = logging.StreamHandler(stream=sys.stdout)
    formatter = jsonlogger.JsonFormatter(LOG_FORMAT, DATETIME_FORMAT)
    log_stdout_handler.setFormatter(formatter)
    logging.basicConfig(handlers=[log_stdout_handler], level=level)
