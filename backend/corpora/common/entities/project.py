import typing

from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink, DbProjectDataset, DbDatasetContributor, DbContributor


class Project(Entity):
    def __init__(self, **kwargs):
        # TODO: list and set fields
        self.id = kwargs['id']
        self.status = kwargs['status']

    @classmethod
    def get(cls, key: typing.Tuple[str]) -> "Project":
        # db_project = cls.db._get(DbProject, key)
        # TODO: query multiple tables: DbProjectLink, DbDataset, etc...
        data = cls.db._query(
            table_args=[DbProject, DbProjectLink],
            filter_args=[
                DbProject.id == key[0],
                DbProject.status == key[1],
                DbProject.id == DbProjectLink.project_id,
                DbProject.status == DbProjectLink.project_status
            ]
        )
        # TODO: inst Project object
        project = Project(id=data.result.DbProject.id,
                          status=data.result.DbProject.status,
                          )
        return project

    @classmethod
    def list(cls):
        pass

    def save(self):
        pass
