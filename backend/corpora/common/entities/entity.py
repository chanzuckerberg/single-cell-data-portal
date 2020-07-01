import typing

from ..corpora_orm import Base
from ..utils.db_utils import DbUtils


class Entity:
    """
    An abstract base class providing an interface to parse application-level objects to and from their
    database counterparts.

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
    def _load(cls, db_result: typing.List["Entity"]) -> "Entity":
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
        return self.db_object.__getattribute__(name)
