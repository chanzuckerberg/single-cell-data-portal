import logging
from contextlib import contextmanager
from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, session as sql_session
from sqlalchemy import event

from .exceptions import CorporaException
from ..corpora_config import CorporaDbConfig
from ..corpora_orm import Base, DbDatasetProcessingStatus

logger = logging.getLogger(__name__)


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


def clone(model: Base, primary_key: dict = None, **kwargs) -> Base:
    """Clone an arbitrary sqlalchemy model object with new primary keys.
    https://stackoverflow.com/questions/28871406/how-to-clone-a-sqlalchemy-db-object-with-new-primary-key

    :param model: The SQLAlchemy model to clone
    :param primary_key: The new primary key values. If None a new primary key is generated
    :param kwargs: Updates the columns in the cloned model.
    :return: a clone of the original model with any kwargs passed in and the new primary key.
    """
    table = model.__table__
    non_pk_columns = [key for key in table.columns.keys() if key not in table.primary_key]
    data = {column: getattr(model, column) for column in non_pk_columns}
    data.update(kwargs)
    if primary_key:
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
