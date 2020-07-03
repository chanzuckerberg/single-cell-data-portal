import typing
import uuid

from backend.corpora.common.utils.exceptions import CorporaException
from ..corpora_orm import Base, DBSessionMaker


class DbUtils:
    def __init__(self):
        self.session = DBSessionMaker().session()
        self.engine = self.session.get_bind()

    def get(self, table: Base, entity_id: typing.Union[str, typing.Tuple[str]]) -> typing.Union[Base, None]:
        """
        Query a table row by its primary key
        :param table: SQLAlchemy Table to query
        :param entity_id: Primary key of desired row
        :return: SQLAlchemy Table object, None if not found
        """
        return self.session.query(table).get(entity_id)

    def query(self, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
        """
        Query the database using the current DB session
        :param table_args: List of SQLAlchemy Tables to query/join
        :param filter_args: List of SQLAlchemy filter conditions
        :return: List of SQLAlchemy query response objects
        """
        return (
            self.session.query(*table_args).filter(*filter_args).all()
            if filter_args
            else self.session.query(*table_args).all()
        )

    retry_limit = 3  # The number of times to attempt generating a uuid

    def generate_id(self, table: Base, *args):
        """
        Generates a UUID to enter into the specified table.

        :param table: The table to generate the uuid  for
        :param arg: Additional primary keys if needed
        :return:
        """
        # Generate the ID
        retry_attempts = 0
        while retry_attempts < self.retry_limit:
            _id = str(uuid.uuid4())
            key = _id if not args else (_id, *args)
            if self.get(table, key):
                retry_attempts += 1
            else:
                return _id
        else:
            raise CorporaException("UUID generation attempt limit exceeded.")
