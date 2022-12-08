import logging
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker

from backend.common.corpora_config import CorporaDbConfig
from backend.common.utils.exceptions import CorporaException

logger = logging.getLogger(__name__)


class DBSessionMaker:
    def __init__(self, database_uri: str = None):
        if not database_uri:
            database_uri = CorporaDbConfig().database_uri
        self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self._session_maker = sessionmaker(bind=self.engine)

    def make_session(self, **kwargs):
        if not self._session_maker:
            self._session_maker = sessionmaker(bind=self.engine)
        return self._session_maker(**kwargs)


# T0D0: remove hard-coded
_db_session_maker = DBSessionMaker(database_uri="postgresql://postgres:secret@localhost")


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
