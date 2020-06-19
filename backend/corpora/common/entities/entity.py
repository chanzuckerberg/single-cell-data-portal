import typing

from ..corpora_orm import Base
from ..utils.db_utils import DbUtils


class Entity:
    """
    An abstract base class providing an interface to parse
    application-level objects to and from their database counterparts.

    Every application-level object must inherit Entity.
    Examples: Project, Dataset
    """
    db = DbUtils()

    @classmethod
    def get(cls, key: typing.Union[str, typing.Tuple[str, str]]):
        """
        Retrieves an entity from the database given its primary key
        :param key: Simple or composite primary key
        :return: Entity
        """
        return cls._load(cls._parse(cls._query(key)))

    @classmethod
    def _query(cls, key: typing.Union[str, typing.Tuple[str, str]]) -> typing.List[Base]:
        """
        Queries the database for required entity data
        :param key: Simple or composite primary key
        :return: list of query result rows
        """
        raise NotImplementedError()

    @classmethod
    def _parse(cls, db_result: typing.List[Base]) -> typing.Dict:
        """
        Parses a database query result into a flattened dict
        :param db_result: list of query result rows
        :return: dict of Entity params
        """
        raise NotImplementedError()

    @classmethod
    def _load(cls, params: typing.Dict) -> "Entity":
        """
        Loads a dict of values into an Entity instance
        :param params: dict of required Entity fields for instantiation
        :return: Entity
        """
        raise NotImplementedError()

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
        :return: saved Entity
        """
        raise NotImplementedError()

