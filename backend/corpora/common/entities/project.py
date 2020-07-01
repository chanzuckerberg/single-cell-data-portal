import typing

from .entity import Entity
from ..corpora_orm import (
    DbProject,
    ProjectStatus,
)
from ..utils.exceptions import CorporaException


class Project(Entity):
    def __init__(self, db_object: DbProject):
        super().__init__(db_object)

    @classmethod
    def _query(cls, key: typing.Tuple[str, str]) -> typing.List[DbProject]:
        """
        Given a key, queries the database for a project and its relevant entities.

        According to the Entity.get interface, the return value of this function is fed as the input to Project._load.

        :param key: Composite primary key tuple of the form (project.id, project.status)
        :return: list of SQLAlchemy query results
        """
        statuses = [status.name for status in ProjectStatus]

        if key[1] not in statuses:
            raise CorporaException(f"Invalid status {key[1]}. Status must be one of {statuses}.")

        result = cls.db.query(table_args=[DbProject], filter_args=[DbProject.id == key[0], DbProject.status == key[1]],)

        return result
