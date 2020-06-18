import typing

from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink, DbProjectDataset, DbDatasetContributor, DbContributor


class Project(Entity):
    def __init__(self, **kwargs):
        # TODO: list and set fields
        self.id = kwargs['id']
        self.status = kwargs['status']
        self.link = kwargs['link']

    @classmethod
    def get(cls, key: typing.Tuple[str, str]) -> "Project":
        # TODO: query multiple tables: DbProjectLink, DbDataset, etc...
        results = cls.db._query(
            table_args=[DbProject, DbProjectLink],
            filter_args=[
                DbProject.id == key[0],
                DbProject.status == key[1],
                DbProject.id == DbProjectLink.project_id,
                DbProject.status == DbProjectLink.project_status
            ]
        )
        # TODO: inst Project object
        project = Project(id=results[0].DbProject.id,
                          status=results[0].DbProject.status,
                          link=results[0].DbProjectLink.link_url
                          )
        return project

    @classmethod
    def list(cls):
        pass

    def save(self):
        pass
