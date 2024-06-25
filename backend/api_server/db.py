from functools import wraps

from flask import g

from backend.common.utils.db_session import db_session_manager

# TODO remove from code base


def dbconnect(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with db_session_manager() as session:
            g.db_session = session
            return func(*args, **kwargs)

    return wrapper
