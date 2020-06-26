import typing

from ..corpora_orm import (
    Base,
    DbProject,
    DbProjectLink,
    DbProjectDataset,
    DbDatasetContributor,
    DbContributor,
    ProjectStatus,
)
from .entity import Entity
from ..utils.exceptions import CorporaException


class Project(Entity):
    def __init__(
        self,
        id: str,
        status: str,
        name: str = "",
        description: str = "",
        owner: str = "",
        s3_bucket: str = "",
        tc_uri: str = "",
        needs_attestation: bool = False,
        processing_state: str = "",
        validation_state: str = "",
        dataset_ids: list = None,
        links: list = None,
        contributors: list = None,
        created_at: str = "",
        updated_at: str = "",
    ):
        self.id = id
        self.status = status
        self.name = name
        self.description = description
        self.owner = owner
        self.s3_bucket = s3_bucket
        self.tc_uri = tc_uri
        self.needs_attestation = needs_attestation
        self.processing_state = processing_state
        self.validation_state = validation_state
        self.dataset_ids = dataset_ids
        self.links = links
        self.contributors = contributors
        self.created_at = created_at
        self.updated_at = updated_at

    @classmethod
    def _query(cls, key: typing.Tuple[str, str]) -> typing.List[Base]:
        """
        Given a key, queries the database for a project and its relevant entities.

        According to the Entity.get interface, the return value of this function is fed as the input to Project._load.

        :param key: Composite primary key tuple of the form (project.id, project.status)
        :return: list of SQLAlchemy query results
        """
        statuses = [status.name for status in ProjectStatus]

        if key[1] not in statuses:
            raise CorporaException(f"Invalid status {key[1]}. Status must be one of {statuses}.")

        result = cls.db.query(
            table_args=[DbProject, DbProjectLink, DbProjectDataset, DbDatasetContributor, DbContributor],
            filter_args=[
                DbProject.id == key[0],
                DbProject.status == key[1],
                DbProject.id == DbProjectLink.project_id,
                DbProject.status == DbProjectLink.project_status,
                DbProject.id == DbProjectDataset.project_id,
                DbProject.status == DbProjectDataset.project_status,
                DbProjectDataset.dataset_id == DbDatasetContributor.dataset_id,
                DbContributor.id == DbDatasetContributor.contributor_id,
            ],
        )

        return result

    @classmethod
    def _load(cls, query_results: typing.List[Base]) -> typing.Union["Project", None]:
        """
        Parses query result rows produced by Project._query into an instantiated Project object.
        The output of this function is the return value of Entity.get.

        SQLAlchemy query results are stored in a list in which each item contains Table objects
        (Db* objects from corpora_orm.py) returned by the query.

        Example query results access:
        row = query_results[0]
        project_id = row.DbProject.id

        :param query_results: list of query result rows
        :return: Project
        """
        if not query_results:
            return None

        # de-dupe peripheral entities from query results
        dataset_ids = set()
        links = {}
        contributors = {}
        for result in query_results:
            dataset_ids.add(result.DbProjectDataset.dataset_id)

            links[result.DbProjectLink.id] = {
                "id": result.DbProjectLink.id,
                "type": result.DbProjectLink.link_type,
                "url": result.DbProjectLink.link_url,
            }

            contributors[result.DbContributor.id] = {
                "id": result.DbContributor.id,
                "name": result.DbContributor.name,
                "institution": result.DbContributor.institution,
                "email": result.DbContributor.email,
            }

        # build dict of Project parameters
        project_params = query_results[0].DbProject.__dict__.copy()
        project_params.pop("_sa_instance_state")
        project_params["dataset_ids"] = list(dataset_ids)
        project_params["links"] = list(links.values())
        project_params["contributors"] = list(contributors.values())

        # instantiate Project
        return Project(**project_params)

    @classmethod
    def list(cls):
        pass

    def save(self):
        pass
