import logging
import typing
from contextlib import contextmanager

from sqlalchemy.exc import SQLAlchemyError

from .exceptions import CorporaException
from ..corpora_orm import Base, DBSessionMaker

logger = logging.getLogger(__name__)

COMMIT_ERROR_MSG = "Failed to commit."


@contextmanager
def session_scope():
    """Database Sessions context manager"""
    session = DBSessionMaker().session()
    try:
        yield session
        logger.info("Commit database changes.")
        session.commit()
    except SQLAlchemyError:
        session.rollback()
        logger.exception(COMMIT_ERROR_MSG)
        raise CorporaException(COMMIT_ERROR_MSG)
    finally:
        logger.info("Closing database session.")
        session.close()


def get(session, table: Base, entity_id: typing.Union[str, typing.Tuple[str]]) -> typing.Union[Base, None]:
    """
    Query a table row by its primary key
    :param session: an open session to the database
    :param table: SQLAlchemy Table to query
    :param entity_id: Primary key of desired row
    :return: SQLAlchemy Table object, None if not found
    """
    return session.query(table).get(entity_id)


def query(session, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
    """
    Query the database using the current DB session
    :param session: an open session to the database
    :param table_args: List of SQLAlchemy Tables to query/join
    :param filter_args: List of SQLAlchemy filter conditions
    :return: List of SQLAlchemy query response objects
    """
    return session.query(*table_args).filter(*filter_args).all() if filter_args else session.query(*table_args).all()


def commit(session):
    """
    Commit changes to the database and roll back if error.

    :param session: an open session to the database
    """
    try:
        session.commit()
    except SQLAlchemyError:
        session.rollback()
        logger.exception(COMMIT_ERROR_MSG)
        raise CorporaException(COMMIT_ERROR_MSG)
