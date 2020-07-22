import logging
import typing
import uuid
from datetime import datetime

from sqlalchemy import and_

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
    list_attributes: typing.Tuple = None  # A list of attributes to retrieve when listing entities
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

    @classmethod
    def list_attributes_in_time_range(
        cls, to_date: int = None, from_date: int = None, filters: list = None, list_attributes: list = None
    ) -> typing.List[typing.Dict]:
        """
        Queries the database for Entities that have been created within the specified time range. Return only the
        entity attributes in `list_attributes`.

        :param to_date: If provided, only lists projects that were created before this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param from_date: If provided, only lists projects that were created after this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param filters: additional filters to apply to the query.
        :param list_attributes: A list of entity attributes to return. If None, the class default is used.
        :return: The results is a list of flattened dictionaries containing the `list_attributes`
        """

        filters = filters if filters else []
        list_attributes = list_attributes if list_attributes else cls.list_attributes
        table = cls.table

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _result[_field] = getattr(db_object, _field)
            return _result

        if to_date:
            filters.append(cls.table.created_at <= datetime.fromtimestamp(to_date))
        if from_date:
            filters.append(table.created_at >= datetime.fromtimestamp(from_date))

        results = [
            to_dict(result)
            for result in cls.db.session.query(table).with_entities(*list_attributes).filter(and_(*filters)).all()
        ]

        return results

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
        Create `rows` in `db_table` associated with Entity Object during object creation. A new UUID is generated and a
        new row is created for each item in `rows`.

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

       :return: a list of database objects to create.
        """
        add_columns = add_columns if add_columns else {}
        db_objs = []
        for columns in rows:
            #  if there are matching keys in columns and add_columns, the key value in add_columns will be used.
            _columns = dict(**columns)
            _columns.update(**add_columns)

            _columns["id"] = str(uuid.uuid4())
            row = db_table(**_columns)
            db_objs.append(row)
        return db_objs
