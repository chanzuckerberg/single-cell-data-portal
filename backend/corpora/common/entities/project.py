import typing

from ..corpora_orm import (
    Base,
    DbProject,
    DbProjectLink,
    DbProjectDataset,
    DbDatasetContributor,
    DbContributor
)
from .entity import Entity


class Project(Entity):
    def __init__(self, **kwargs):
        self.id = kwargs['id']
        self.status = kwargs['status']
        self.name = kwargs.get("name", "")
        self.description = kwargs.get("description", "")
        self.submitter = kwargs.get('submitter', "")
        self.s3_bucket = kwargs.get('s3_bucket', "")
        self.tc_uri = kwargs.get('tc_uri', "")
        self.needs_attestation = kwargs.get('needs_attestation', "")
        self.processing_state = kwargs.get('processing_state', "")
        self.validation_state = kwargs.get('validation_state', "")
        self.dataset_ids = kwargs.get('dataset_ids', [])
        self.links = kwargs.get('links', [])
        self.contributors = kwargs.get('contributors', [])

    @classmethod
    def _query(cls, key: typing.Tuple[str, str]) -> typing.List[Base]:
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
                DbContributor.id == DbDatasetContributor.contributor_id
            ]
        )
        return result

    @classmethod
    def _parse(cls, query_results):
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

    @classmethod
    def _load(cls, params: typing.Dict):
        return Project(id=params['id'],
                       status=params['status'],
                       name=params['name'],
                       description=params['description'],
                       submitter=params['owner'],
                       s3_bucket=params['s3_bucket'],
                       tc_uri=params['tc_uri'],
                       needs_attestation=params['needs_attestation'],
                       processing_state=params['processing_state'],
                       validation_state=params['validation_state'],
                       dataset_ids=params['dataset_ids'],
                       links=params['links'],
                       contributors=params['contributors'])

    @classmethod
    def list(cls):
        pass

    def save(self):
        pass
