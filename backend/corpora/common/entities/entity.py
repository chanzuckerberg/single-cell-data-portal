import logging
import typing
import uuid

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
    def list(cls) -> "Entity":
        """
        Retrieves a list of entities from the database
        :return: list of Entity
        """
        return [cls(obj) for obj in cls.db.query([cls.table])]

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

    @classmethod
    def _create_sub_objects(
        cls, rows: typing.List[dict], db_table: Base, add_columns: dict = None, primary_keys: typing.List[str] = None
    ) -> typing.List[Base]:
        """
        Create or modify `rows` in `db_table` associated with Entity Object during object creation. If id is provided,
        then the row with that id is retrieved and updated, else a new UUID is generated and a new row is created.

        :param rows: A list of dictionaries each specifying a row to insert or modify
        :param db_table: The Table to add or modify rows
        :param primary_keys: Additional columns required to build the primary key. This is used when the primary consist
        of multiple columns.
        :param add_columns: Additional columns attributes or modifications to add to the row.

        This can be used when there are shared column values that need to be added across all the new rows.
        For example: DbProjectLink generated for a specific project should all have the same DbProjectLink.project_id
        and DbProjectLink.project_status. The function call would be:
        >>>> cls._create_sub_objects(
        >>>>    [{'link_url':'abc', 'link_type': ProjectLinkType.OTHER}],
        >>>>    DbProjectLink,
        >>>>    add_columns={'project_id':'abcd','project_status':ProjectStatus.EDIT}
        >>>>    )

        Another use would be to overwrite column specified in the rows.

       :return: a list of database objects to create or modify.
        """
        add_columns = add_columns if add_columns else {}
        db_objs = []
        for columns in rows:
            _columns = dict(**columns)
            _columns.update(**add_columns)  # if there are matching keys in columns and add_columns,
            # the key value in add_columns will be used.
            _id = _columns.get("id")
            if _id:
                primary_key = (_id, *primary_keys) if primary_keys else _id
                row = cls.db.get(db_table, primary_key)
            else:
                _columns["id"] = str(uuid.uuid4())
                row = db_table(**_columns)
            db_objs.append(row)
        return db_objs
