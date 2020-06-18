import typing

from ..corpora_orm import Base, DBSessionMaker, DbProject


class DbUtils:
    def __init__(self):
        self.session = DBSessionMaker().session()
        self.engine = self.session.get_bind()

    def _get(self, table: Base, entity_id: typing.Union[str, typing.Tuple[str]]) -> typing.Union[Base, None]:
        """
        Query a table row by its primary key
        :param table: SQLAlchemy Table to query
        :param entity_id: Primary key of desired row
        :return: SQLAlchemy Table object, None if not found
        """
        return self.session.query(table).get(entity_id)

    def _query(self, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
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

    @staticmethod
    def _parse_multivalue(value: str) -> typing.List[str]:
        """
        Parses a CSV string representing multiple values into a list
        :param value: Comma-separated value to parse
        :return: List of strings
        """
        return value.split(",") if value else []

    def query_projects(self):
        return [
            {
                "id": project.id,
                "status": project.status,
                "name": project.name,
                "description": project.description,
                "processing_state": project.processing_state.name,
                "validation_state": project.validation_state.name,
            }
            for project in self._query([DbProject])
        ]
