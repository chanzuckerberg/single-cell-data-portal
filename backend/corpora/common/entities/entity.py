import logging
import typing

logger = logging.getLogger(__name__)

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

    table: Base = None  # The DbTable represented by this entity.
    db = DbUtils()

    def __init__(self, db_object: Base):
        self.db_object = db_object

    @classmethod
    def get(cls, key: typing.Union[str, typing.Tuple[str, str]]) -> typing.Union["Entity", None]:
        """
        Retrieves an entity from the database given its primary key if found.
        :param key: Simple or composite primary key
        :return: Entity or None
        """
        result = cls.db.get(cls.table, key)
        if result:
            return cls(result)
        else:
            logger.info(f"Unable to find a row with primary key {key}, in {cls.__name__} table.")
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

    @classmethod
    def _create_sub_objects(
        cls, rows: typing.List[dict], db_table: Base, keys: typing.List[str] = None, **kwargs
    ) -> typing.List[Base]:
        """
        Create a list of Table Rows to be added to an Entity Object during an object creation. If id is provided, then
        the row is retrieved and updated.

        :param rows: A list of dictionaries containing columns.
        :param db_table: The Table to add or modify rows
        :param keys: additional primary keys
        :param kwargs: Additional columns attributes to add.
        :return:
        """
        db_objs = []
        for columns in rows:
            _columns = dict(**columns, **kwargs)
            _id = _columns.get("id")
            if _id:
                primary_key = (_id, *keys) if keys else _id
                row = cls.db.get(db_table, primary_key)
            else:
                _columns["id"] = cls.db.generate_id(db_table)
                row = db_table(**_columns)
            db_objs.append(row)
        return db_objs
