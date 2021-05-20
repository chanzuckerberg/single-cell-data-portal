import logging
from contextlib import contextmanager
from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, session

from .exceptions import CorporaException
from ..corpora_config import CorporaDbConfig

logger = logging.getLogger(__name__)


class DBSessionMaker:

    _session_make = None
    engine = None

    def __init__(self, database_uri: str = None):
        if not self.engine:
            database_uri = database_uri if database_uri else CorporaDbConfig().database_uri
            self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        if not self._session_make:
            self._session_make = sessionmaker(bind=self.engine)

    def session(self, **kwargs) -> session.Session:
        return self._session_make(**kwargs)


@contextmanager
def db_session_manager(**kwargs):
    """

    :param kwargs: passed to Session
    """
    try:
        session = DBSessionMaker().session(**kwargs)
        yield session
        if session.transaction:
            session.commit()
        else:
            session.expire_all()
    except SQLAlchemyError:
        session.rollback()
        msg = "Failed to commit."
        logger.exception(msg)
        raise CorporaException(msg)
    finally:
        session.close()
