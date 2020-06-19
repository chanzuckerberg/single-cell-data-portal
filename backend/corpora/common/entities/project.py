import typing

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
    def _query(cls, key: typing.Tuple[str, str]):
        return cls.db.get_project(key)

    @classmethod
    def _parse(cls, query_results):
        return cls.db.parse_project(query_results)

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
