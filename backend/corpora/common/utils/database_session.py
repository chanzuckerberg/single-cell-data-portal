import functools
import logging
from contextlib import contextmanager

import sqlalchemy
from sqlalchemy.exc import SQLAlchemyError

from .exceptions import CorporaException
from .singleton import Singleton
from ..corpora_orm import DBSessionMaker

logger = logging.getLogger(__name__)


class DatabaseSession(metaclass=Singleton):
    """DbUtils as a singleton to avoid creating excess sessions."""

    def __init__(self):
        self._session = None
        self.session_maker = DBSessionMaker().session_maker

    @property
    def session(self) -> sqlalchemy.orm.session.Session:
        if not self._session:
            self._session = self.session_maker()
        return self._session

    def __getattr__(self, item):
        return getattr(self.session, item)


@contextmanager
def db_session_manager(commit=False):
    """

    :param commit: Changes will be committed when context ends.
    """
    try:
        session = DatabaseSession()
        yield session
        if commit:
            session.commit()
    except SQLAlchemyError:
        session.rollback()
        msg = "Failed to commit."
        logger.exception(msg)
        raise CorporaException(msg)
    finally:
        session.close()


def db_session(commit=False):
    """

    :param commit: passed to db_session_manager
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with db_session_manager(commit):
                rv = func(*args, **kwargs)
            return rv

        return wrapper

    return decorator
