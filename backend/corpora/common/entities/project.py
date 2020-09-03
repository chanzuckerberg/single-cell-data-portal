import typing
import uuid
from datetime import datetime

from sqlalchemy import and_

from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink, ProjectStatus


class Project(Entity):
    table = DbProject
    list_attributes = (DbProject.id, DbProject.created_at)

    def __init__(self, db_object: DbProject):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        status: str,
        name: str = "",
        description: str = "",
        owner: str = "",
        s3_bucket: str = "",
        tc_uri: str = "",
        needs_attestation: bool = False,
        processing_state: str = "",
        validation_state: str = "",
        links: list = None,
        **kwargs,
    ) -> "Project":
        """
        Create a new Project and related objects and store in the database. UUIDs are generated for all new table
        entries.
        """
        primary_key = str(uuid.uuid4())

        # Setting Defaults
        links = links if links else []

        new_db_object = DbProject(
            id=primary_key,
            status=status,
            name=name,
            description=description,
            owner=owner,
            s3_bucket=s3_bucket,
            tc_uri=tc_uri,
            needs_attestation=needs_attestation,
            processing_state=processing_state,
            validation_state=validation_state,
            links=cls._create_sub_objects(
                links, DbProjectLink, add_columns=dict(project_id=primary_key, project_status=status)
            ),
            **kwargs,
        )

        cls.db.session.add(new_db_object)
        cls.db.commit()
        return cls(new_db_object)

    @classmethod
    def get_project(cls, project_uuid):
        """
        Given the project_uuid, retrieve a live project.
        :param project_uuid:
        """
        return cls.get((project_uuid, ProjectStatus.LIVE.name))

    @classmethod
    def list_projects_in_time_range(cls, *args, **kwargs):
        return cls.list_attributes_in_time_range(*args, filters=[DbProject.status == ProjectStatus.LIVE.name], **kwargs)

    @classmethod
    def list_attributes_in_time_range(
        cls, to_date: int = None, from_date: int = None, filters: list = None, list_attributes: list = None
    ) -> typing.List[typing.Dict]:
        """
        Queries the database for Entities that have been created within the specified time range. Return only the
        entity attributes in `list_attributes`.

        :param to_date: If provided, only lists projects that were created before this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param from_date: If provided, only lists projects that were created after this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param filters: additional filters to apply to the query.
        :param list_attributes: A list of entity attributes to return. If None, the class default is used.
        :return: The results is a list of flattened dictionaries containing the `list_attributes`
        """

        filters = filters if filters else []
        list_attributes = list_attributes if list_attributes else cls.list_attributes
        table = cls.table

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _result[_field] = getattr(db_object, _field)
            return _result

        if to_date:
            filters.append(cls.table.created_at <= datetime.fromtimestamp(to_date))
        if from_date:
            filters.append(table.created_at >= datetime.fromtimestamp(from_date))

        results = [
            to_dict(result)
            for result in cls.db.session.query(table).with_entities(*list_attributes).filter(and_(*filters)).all()
        ]

        return results

    @classmethod
    def get_submission(cls, project_uuid):
        """
        Given the project_uuid, retrieve a live project.
        :param project_uuid:
        """
        return cls.get((project_uuid, ProjectStatus.EDIT.name))

    @classmethod
    def list_submissions(cls, *args, **kwargs):
        return cls.list_attributes_in_time_range(
            *args,
            filters=[DbProject.status == ProjectStatus.EDIT.name],
            list_attributes=[
                DbProject.id,
                DbProject.name,
                DbProject.processing_state,
                DbProject.validation_state,
                DbProject.owner,
            ],
            **kwargs,
        )

    def reshape_for_api(self) -> dict:
        """
        Reshape the project to match the expected api output.
        :return: A dictionary that can be converted into JSON matching the expected api response.
        """
        result = self.to_dict()
        # Reshape the data to match.
        result["s3_bucket_key"] = result.pop("s3_bucket", None)
        result.pop("user", None)
        result.pop("owner", None)
        result["links"] = [
            dict(url=link["link_url"], name=link["link_name"], type=link["link_type"]) for link in result["links"]
        ]
        result["attestation"] = dict(needed=result.pop("needs_attestation", None), tc_uri=result.pop("tc_uri", None))
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = dataset.pop("deployment_directories")
            dataset["dataset_assets"] = dataset.pop("artifacts")
            dataset.pop("contributors", None)
        return result
