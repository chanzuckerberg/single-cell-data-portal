import typing
import uuid
from datetime import datetime

import pytz
from sqlalchemy import and_

from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink, ProjectStatus


class Project(Entity):
    table = DbProject

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

        #  Prevent accidentally linking an existing row to a different Project. This maintains the relationship of one
        #  to many for links
        [link.pop("id", None) for link in links]

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
                links, DbProjectLink, add_columns=dict(project_id=primary_key, project_status=status, **kwargs)
            ),
            **kwargs
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
        return cls.list_in_time_range(*args, filters=[DbProject.status == ProjectStatus.LIVE.name], **kwargs)

    @classmethod
    def list_in_time_range(
        cls, to_date: float = None, from_date: float = None, filters: list = None
    ) -> typing.List[typing.Dict]:
        """
        Queries the database for projects that have been created within the specified time range.

        :param to_date: If provided, only lists projects that were created before this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param from_date: If provided, only lists projects that were created after this date. Format of param is Unix timestamp since the epoch in UTC timezone.
        :param filters: additional filters to apply to the query.
        :return: The results is a list of flattened dictionaries containing the `list_entities`
        """

        filters = filters if filters else []
        list_entities = [DbProject.created_at, DbProject.id]

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _result[_field] = getattr(db_object, _field)
            return _result

        if to_date:
            filters.append(DbProject.created_at <= datetime.fromtimestamp(to_date, tz=pytz.UTC))
        if from_date:
            filters.append(DbProject.created_at >= datetime.fromtimestamp(from_date, tz=pytz.UTC))

        results = [
            to_dict(result)
            for result in cls.db.session.query(DbProject).with_entities(*list_entities).filter(and_(*filters)).all()
        ]

        return results
