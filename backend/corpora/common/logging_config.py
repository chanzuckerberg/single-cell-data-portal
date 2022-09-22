LOGGED_FIELDS = ["levelname", "asctime", "name", "message", "lineno", "pathname"]
LOG_FORMAT = " ".join([f"%({field})" for field in LOGGED_FIELDS])

DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%03dZ"
