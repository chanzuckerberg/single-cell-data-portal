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
    def __init__(self, **kwargs):
        self.id = kwargs["id"]
        self.status = kwargs["status"]
        self.name = kwargs.get("name", "")
        self.description = kwargs.get("description", "")
        self.submitter = kwargs.get("submitter", "")
        self.s3_bucket = kwargs.get("s3_bucket", "")
        self.tc_uri = kwargs.get("tc_uri", "")
        self.needs_attestation = kwargs.get("needs_attestation", "")
        self.processing_state = kwargs.get("processing_state", "")
        self.validation_state = kwargs.get("validation_state", "")
        self.dataset_ids = kwargs.get("dataset_ids", [])
        self.links = kwargs.get("links", [])
        self.contributors = kwargs.get("contributors", [])

    @classmethod
    def _query(cls, key: typing.Tuple[str, str]) -> typing.List[Base]:
        """
        Given a key, queries the database for a project and its relevant entities.

        According to the Entity.get interface, the return value of this function is
        fed as the input to Project._parse.

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
    def _parse(cls, query_results: typing.List[Base]) -> typing.Dict:
        """
        Parses query result rows produced by Project._query into a dict of KVPs
        required for entity instantiation.
        The output of this function is the input to Project._load.

        SQLAlchemy query results are stored in a list in which each item
        contains Table objects (Db* objects from corpora_orm.py) returned by the query.

        Example query results access:
        row = query_results[0]
        project_id = row.DbProject.id

        :param query_results: list of query result rows
        :return: dict of KVPs
        """
        if not query_results:
            return {}

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
        project = {}
        for k, v in query_results[0].DbProject.__dict__.items():
            project[k] = v
        project["dataset_ids"] = list(dataset_ids)
        project["links"] = list(links.values())
        project["contributors"] = list(contributors.values())

        return project

    @classmethod
    def _load(cls, params: typing.Dict):
        """
        Instantiates and returns a Project instance.

        :param params: dict of KVPs used to instantiate a Project
        :return: Project
        """
        return Project(
            id=params["id"],
            status=params["status"],
            name=params["name"],
            description=params["description"],
            submitter=params["owner"],
            s3_bucket=params["s3_bucket"],
            tc_uri=params["tc_uri"],
            needs_attestation=params["needs_attestation"],
            processing_state=params["processing_state"],
            validation_state=params["validation_state"],
            dataset_ids=params["dataset_ids"],
            links=params["links"],
            contributors=params["contributors"],
        ) if params else None

    @classmethod
    def list(cls):
        pass

    def save(self):
        pass
