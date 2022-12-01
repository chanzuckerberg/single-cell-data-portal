import logging

from contextlib import contextmanager
from flask import g
from functools import wraps

from sqlalchemy import create_engine
from sqlalchemy import event
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, session as sql_session

from backend.common.utils.exceptions import CorporaException
from backend.common.corpora_config import CorporaDbConfig

logger = logging.getLogger(__name__)


class DBSessionMaker:

    def __init__(self, database_uri: str = None):
        if not database_uri:
            database_uri = CorporaDbConfig().database_uri
            # database_uri = "postgresql://postgres:secret@localhost"
        self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self._session_maker = sessionmaker(bind=self.engine)

    def make_session(self, **kwargs):
        if not self._session_maker:
            self._session_maker = sessionmaker(bind=self.engine)
        return self._session_maker(**kwargs)


_db_session_maker = DBSessionMaker()


@contextmanager
def db_session_manager(db_session_maker=_db_session_maker, **kwargs):
    try:
        session = db_session_maker.make_session(**kwargs)
        yield session
        if session.transaction:
            session.commit()
        else:
            session.expire_all()
    except SQLAlchemyError as e:
        logger.exception(e)
        if session is not None:
            session.rollback()
        raise CorporaException("Failed to commit.")
    finally:
        if session is not None:
            session.close()
