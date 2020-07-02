import typing

from ..corpora_orm import Base
from ..utils.db_utils import DbUtils


class Entity:
    """
    An abstract base class providing an interface to parse application-level objects to and from their
    database counterparts.

    This class uses a has-a relationship with SQLAlchemy Table object and simplify the CRUD operations performed on the
    database through these objects. Columns and relationships of the database object can be access as attributes of the
    of Entity.

    Every application-level object must inherit Entity.
    Examples: Project, Dataset
    """

    db = DbUtils()

    def __init__(self, db_object: Base):
        self.db_object = db_object

    @classmethod
    def get(cls, key: typing.Union[str, typing.Tuple[str, str]]):
        """
        Retrieves an entity from the database given its primary key
        :param key: Simple or composite primary key
        :return: Entity
        """
        return cls._load(cls._query(key))

    @classmethod
    def _query(cls, key: typing.Union[str, typing.Tuple[str, str]]) -> typing.List[Base]:
        """
        Queries the database for required entity data
        :param key: Simple or composite primary key
        :return: list of query result rows
        """
        raise NotImplementedError()

    @classmethod
    def _load(cls, db_result: typing.List["Entity"]) -> typing.Union["Entity", None]:
        """

        Parses query result rows produced by Entity._query into an instantiated Entity object.
        The output of this function is the return value of Entity.get.
        SQLAlchemy query results are stored in a list in which each item contains Table objects returned by the query.
        If no results are found, None is returned.

        Example query results access:
        row = query_results[0]
        project_id = row.DbProject.id
        :param query_results: list of query result rows
        :return: Entity
        """

        try:
            return cls(db_result[0])
        except IndexError:
            return None

    @classmethod
    def list(cls):
        """
        Retrieves a list of entities from the database
        :return: list of Entity
        """
        raise NotImplementedError()

    def save(self):
        """
        Writes the current object state to the database
        :return: saved Entity object
        """
        raise NotImplementedError()

    def __getattr__(self, name):
        """
        If the attribute is not in Entity, return the attribute in db_object.
        :param name:
        :return:
        """
        return self.db_object.__getattribute__(name)
