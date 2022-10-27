from typing import List


def format_log_string(fields: List[str]) -> str:
    return " ".join([f"%({field})" for field in fields])


LOGGED_FIELDS = ["levelname", "asctime", "name", "message", "lineno", "pathname"]
LOG_FORMAT = format_log_string(LOGGED_FIELDS)


DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%03dZ"
