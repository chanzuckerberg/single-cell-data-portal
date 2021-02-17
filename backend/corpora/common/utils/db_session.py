import logging
from contextlib import contextmanager
from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, session

from .exceptions import CorporaException
from ..corpora_config import CorporaDbConfig
from ..corpora_orm import Base, DbDatasetProcessingStatus

logger = logging.getLogger(__name__)


class DBSessionMaker:

    _session_make = None
    engine = None

    def __init__(self):
        if not self.engine:
            self.engine = create_engine(CorporaDbConfig().database_uri, connect_args={"connect_timeout": 5})
        if not self._session_make:
            self._session_make = sessionmaker(bind=self.engine)

    def session(self, **kwargs) -> session.Session:
        return self._session_make(**kwargs)


def clone(model: Base, primary_key: dict, **kwargs) -> Base:
    """Clone an arbitrary sqlalchemy model object with new primary keys.
    https://stackoverflow.com/questions/28871406/how-to-clone-a-sqlalchemy-db-object-with-new-primary-key

    :param model: The SQLAlchemy model to clone
    :param primary_key: The new primary key values.
    :param kwargs: Updates the columns in the cloned model.
    :return: a clone of the original model with any kwargs passed in and the new primary key.
    """
    table = model.__table__
    non_pk_columns = [key for key in table.columns.keys() if key not in table.primary_key]
    data = {column: getattr(model, column) for column in non_pk_columns}
    data.update(kwargs)
    data.update(primary_key)
    return model.__class__(**data)


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


def processing_status_updater(session, uuid: str, updates: dict):
    session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == uuid).update(updates)
    session.commit()
    logger.debug("updating status", updates)
