import logging
import sys

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT


def configure_logging():
    log_stdout_handler = logging.StreamHandler(stream=sys.stdout)
    formatter = jsonlogger.JsonFormatter(LOG_FORMAT, DATETIME_FORMAT)
    log_stdout_handler.setFormatter(formatter)
    logging.basicConfig(handlers=[log_stdout_handler], level=logging.INFO)


def logit(func):
    def wrapper(*arg, **kw):
        """Logging the start and finish of a function"""
        logging.info(f"Start {func.__name__}", extra={"type": "METRIC"})
        res = func(*arg, **kw)
        logging.info(f"Complete {func.__name__}", extra={"type": "METRIC"})
        return res

    return wrapper
