import typing

from ..corpora_orm import (
    Base,
    DBSessionMaker,
    DbProject,
    DbProjectLink,
    DbProjectDataset,
    DbDatasetContributor,
    DbContributor
)


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

    def get_project(self, key: typing.Tuple[str, str]) -> typing.List[Base]:
        result = self._query(
            table_args=[DbProject, DbProjectLink, DbProjectDataset, DbDatasetContributor, DbContributor],
            filter_args=[
                DbProject.id == key[0],
                DbProject.status == key[1],
                DbProject.id == DbProjectLink.project_id,
                DbProject.status == DbProjectLink.project_status,
                DbProject.id == DbProjectDataset.project_id,
                DbProject.status == DbProjectDataset.project_status,
                DbProjectDataset.dataset_id == DbDatasetContributor.dataset_id,
                DbContributor.id == DbDatasetContributor.contributor_id
            ]
        )
        return result

    @staticmethod
    def parse_project(query_results: typing.List[Base]):
        # de-dupe and parse project relationships
        dataset_ids = set()
        links = {}
        contributors = {}
        for result in query_results:
            dataset_ids.add(result.DbProjectDataset.dataset_id)

            links[result.DbProjectLink.id] = {
                'id': result.DbProjectLink.id,
                'type': result.DbProjectLink.link_type,
                'url': result.DbProjectLink.link_url
            }

            contributors[result.DbContributor.id] = {
                'id': result.DbContributor.id,
                'name': result.DbContributor.name,
                'institution': result.DbContributor.institution,
                'email': result.DbContributor.email
            }

        # build project params
        project = {}
        for k, v in query_results[0].DbProject.__dict__.items():
            project[k] = v
        project['dataset_ids'] = list(dataset_ids)
        project['links'] = list(links.values())
        project['contributors'] = list(contributors.values())

        return project
