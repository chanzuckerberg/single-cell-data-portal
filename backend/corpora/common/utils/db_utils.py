import functools
import logging
import typing
from contextlib import contextmanager

import sqlalchemy
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

        def get(self, table: Base, entity_id: typing.Union[str, typing.Tuple[str]]) -> typing.Union[Base, None]:
            """
            Query a table row by its primary key
            :param table: SQLAlchemy Table to query
            :param entity_id: Primary key of desired row
            :return: SQLAlchemy Table object, None if not found
            """
            return self.session.query(table).get(entity_id)

        def query(self, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
            """
            Query the database using the current DB session
            :param table_args: List of SQLAlchemy Tables to query/join
            :param filter_args: List of SQLAlchemy filter conditions
            :return: List of SQLAlchemy query response objects
            """
            return (
                self.session.query(*table_args).filter(*filter_args).all()
                if filter_args
                else self.session.query(*table_args).all()
            )

        def commit(self):
            """
            Commit changes to the database and roll back if error.
            """
            self.session.commit()

        def delete(self, db_object: Base):
            try:
                self.session.delete(db_object)
            except sqlalchemy.orm.exc.ObjectDeletedError:
                pass

        def close(self):
            self.session.close()
            self._session = None

    def __init__(self):
        if not DbUtils.__instance:
            DbUtils.__instance = DbUtils.__DbUtils()

    def __getattr__(self, name):
        return getattr(self.__instance, name)


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
def db_session_manager(commit=False):
    """

    :param commit: Changes will be committed when context ends.
    """
    try:
        db = DbUtils()
        yield db
        if commit:
            db.commit()
    except SQLAlchemyError:
        db.session.rollback()
        msg = "Failed to commit."
        logger.exception(msg)
        raise CorporaException(msg)
    finally:
        db.close()


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


def processing_status_updater(uuid: str, updates: dict):
    with db_session_manager(commit=True) as manager:
        manager.session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == uuid).update(updates)
