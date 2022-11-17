from ddtrace import tracer
from flask import g
from functools import wraps

from backend.common.utils.db_session import db_session_manager


@tracer.wrap()
def dbconnect_and_ddtrace(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with db_session_manager() as session:
            g.db_session = session
            return func(*args, **kwargs)

    return wrapper
