import logging
import typing
from datetime import datetime
from operator import and_

from backend.corpora.common.corpora_orm import Base

logger = logging.getLogger(__name__)


def get(session, table, key: typing.Union[str, typing.Tuple[str, str]]) -> typing.Union[Base, None]:
    """
    Retrieves an entity from the database given its primary key if found.
    :param session:
    :param table:
    :param key:
    :return:
    """
    result = session.get(table, key)
    if result:
        return result
    else:
        logger.info(f"Unable to find a row with primary key {key}, in {table} table.")
        return None


def list_table(session, table) -> typing.List[Base]:
    """
    Retrieves a list of entities from the database
    :return: list of Entity
    """
    return [obj for obj in session.query([table])]


def delete(session, db_object):
    """
    Delete an object from the database.
    """
    session.delete(db_object)
    session.commit()


def update(session, db_object, **kwargs):
    for key, value in kwargs.items():
        if hasattr(db_object, key):
            setattr(db_object, key, value)
    session.commit()


def list_attributes_in_time_range(
    session, table, to_date: int = None, from_date: int = None, filters: list = None, list_attributes: list = None
) -> typing.List[typing.Dict]:
    """
    Queries the database for rows that have been created within the specified time range. Return only the
    row attributes in `list_attributes`.

    :param session:
    :param table:
    :param to_date: If provided, only lists collections that were created before this date. Format of param is Unix
    timestamp since the epoch in UTC timezone.
    :param from_date: If provided, only lists collections that were created after this date. Format of param is Unix
    timestamp since the epoch in UTC timezone.
    :param filters: additional filters to apply to the query.
    :param list_attributes: A list of entity attributes to return. If None, the class default is used.
    :return: The results is a list of flattened dictionaries containing the `list_attributes`
    """

    filters = filters if filters else []
    list_attributes = list_attributes if list_attributes else list_attributes

    def to_dict(db_object):
        _result = {}
        for _field in db_object._fields:
            _result[_field] = getattr(db_object, _field)
        return _result

    if to_date:
        filters.append(table.created_at <= datetime.fromtimestamp(to_date))
    if from_date:
        filters.append(table.created_at >= datetime.fromtimestamp(from_date))

    if list_attributes:
        results = [
            to_dict(result)
            for result in session.query(table).with_entities(*list_attributes).filter(and_(*filters)).all()
        ]
    else:
        results = [to_dict(result) for result in session.query(table).filter(and_(*filters)).all()]

    return results
