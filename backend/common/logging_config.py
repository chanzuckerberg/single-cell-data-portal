import logging
from contextlib import suppress
from typing import List


def format_log_string(fields: List[str]) -> str:
    return " ".join([f"%({field})" for field in fields])


LOGGED_FIELDS = ["levelname", "asctime", "name", "message", "lineno", "pathname"]
LOG_FORMAT = format_log_string(LOGGED_FIELDS)


DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%03dZ"


def logit(func):
    def wrapper(*arg, **kw):
        """Logging the start and finish of a function"""
        logging.info(f"Start {func.__name__}", extra={"type": "METRIC"})
        res = func(*arg, **kw)
        logging.info(f"Complete {func.__name__}", extra={"type": "METRIC"})
        return res

    return wrapper


class LogSuppressed(suppress):
    def __init__(self, *args, message="Suppressed Exception"):
        super().__init__(*args)
        self.message = message

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            logging.error(self.message, exc_info=(exc_type, exc_value, traceback))
        return super().__exit__(exc_type, exc_value, traceback)
