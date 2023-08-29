import logging
from contextlib import contextmanager

from sqlalchemy import create_engine, event
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import session as sql_session
from sqlalchemy.orm import sessionmaker

from backend.common.corpora_config import CorporaDbConfig
from backend.common.utils.exceptions import CorporaException

logger: logging.Logger = logging.getLogger(__name__)


class DBSessionMaker:

    _session_maker = None
    engine = None

    def __init__(self, database_uri: str = None):
        if not self.engine:
            self.database_uri = database_uri if database_uri else CorporaDbConfig().database_uri
            self.engine = create_engine(self.database_uri, connect_args={"connect_timeout": 5})

    def session(self, **kwargs) -> sql_session.Session:
        new_session = self.session_maker(info=dict(s3_deletion_list=[]), **kwargs)

        @event.listens_for(new_session, "after_commit")
        def cleanup_s3_objects(session):
            for func in session.info["s3_deletion_list"]:
                func()
            session.info["s3_deletion_list"].clear()

        return new_session

    @property
    def session_maker(self):
        if not self._session_maker:
            self._session_maker = sessionmaker(bind=self.engine)
        return self._session_maker


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
    except SQLAlchemyError as e:
        logger.exception(e)
        if session is not None:
            session.rollback()
        raise CorporaException("Failed to commit.") from None
    finally:
        if session is not None:
            session.close()
