import typing
from abc import ABC, abstractmethod

from ..utils.db_utils import DbUtils


class Entity(ABC):
    """
    An abstract base class providing an interface to parse
    application-level objects to and from their database counterparts.

    Every application-level object must inherit Entity.
    Examples: Project, Dataset
    """
    # Proposal: implement Singleton via ABC class variable to minimize # sessions created.
    # TODO: verify attributes can be set in an ABC
    db = DbUtils()

    @classmethod
    @abstractmethod
    def get(cls, key: typing.Union[str, typing.Tuple[str, str]]):
        """
        Retrieves an entity from the database given its primary key
        :param key: Simple or composite primary key
        :return: Entity
        """
        pass

    @classmethod
    @abstractmethod
    def list(cls):
        """
        Retrieves a list of entities from the database
        :return: list of Entity
        """
        pass

    @abstractmethod
    def save(self):
        """
        Writes the current object state to the database
        :return: saved Entity
        """
        pass
