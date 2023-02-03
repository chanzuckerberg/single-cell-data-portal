import logging
import typing

from sqlalchemy import inspect
from sqlalchemy.orm import Session

from backend.common.corpora_orm import Base


class Entity:
    """
    An abstract base class providing an interface to parse application-level objects to and from their
    database counterparts.

    This class uses a has-a relationship with SQLAlchemy Table object and simplify the CRUD operations performed on the
    database through these objects. Columns and relationships of the database object can be access as attributes of the
    of Entity.

    Every application-level object must inherit Entity.
    Examples: Collection, Dataset
    """

    table: Base = None  # The DbTable represented by this entity.
    list_attributes: typing.Tuple = None  # A list of attributes to retrieve when listing entities

    def __init__(self, db_object: Base):
        self.db_object = db_object
        self.session = inspect(db_object).session

    @classmethod
    def get(cls, session: Session, key: typing.Union[str, typing.Tuple[str, str]]) -> typing.Union["Entity", None]:
        """
        Retrieves an entity from the database given its primary key if found.
        :param session: The SQLAlchemy database Session
        :param key: Simple or composite primary key
        :return: Entity or None
        """
        result = session.query(cls.table).get(key)
        if result:
            return cls(result)
        else:
            logging.info(f"Unable to find a row with primary key {key}, in {cls.__name__} table.")
            return None

    @classmethod
    def list(cls, session: Session) -> "Entity":
        """
        Retrieves a list of entities from the database
        :return: list of Entity
        """
        return [cls(obj) for obj in session.query(cls.table)]

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
        """
        return self.db_object.__getattribute__(name)

    def delete(self, commit=True):
        """
        Delete an object from the database.
        """
        self.session.delete(self.db_object)
        if commit:
            self.session.commit()

    def update(self, commit=True, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self.db_object, key):
                setattr(self.db_object, key, value)
        if commit:
            self.session.commit()

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id})>"
