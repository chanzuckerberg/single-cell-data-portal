import logging
from contextlib import contextmanager


from sqlalchemy.exc import SQLAlchemyError

from .exceptions import CorporaException
from ..corpora_orm import Base, DBSessionMaker, DbDatasetProcessingStatus

logger = logging.getLogger(__name__)


class DbUtils:
    """DbUtils as a singleton to avoid creating excess sessions."""

    __instance = None

    class __DbUtils:
        def __init__(self):
            self._session = None

        @property
        def session(self):
            if not self._session:
                self._session = DBSessionMaker().session()
            return self._session

    def __init__(self):
        if not DbUtils.__instance:
            DbUtils.__instance = DbUtils.__DbUtils()

    def __getattr__(self, name):
        return getattr(self.__instance, name, None) or getattr(self.__instance.session, name)


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
def db_session_manager(commit=False, **kwargs):
    """

    :param commit: Changes will be committed when context ends.
    """
    try:
        session = DBSessionMaker().session(**kwargs)
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


def processing_status_updater(session, uuid: str, updates: dict):
    session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == uuid).update(updates)
    session.commit()
    logger.debug("updating status", updates)
