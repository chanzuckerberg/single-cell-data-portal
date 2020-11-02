import functools
import logging
import typing

from sqlalchemy.exc import SQLAlchemyError

from .exceptions import CorporaException
from ..corpora_orm import Base, DBSessionMaker

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
            self.session.delete(db_object)

        def close(self):
            self.session.close()
            self._session = None

    def __init__(self):
        if not DbUtils.__instance:
            DbUtils.__instance = DbUtils.__DbUtils()

    def __getattr__(self, name):
        return getattr(self.__instance, name)


def db_session(func):
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        db = DbUtils()
        try:
            rv = func(*args, **kwargs)
            return rv
        except SQLAlchemyError:
            db.session.rollback()
            msg = "Failed to commit."
            logger.exception(msg)
            raise CorporaException(msg)
        finally:
            db.close()

    return wrapper_decorator
