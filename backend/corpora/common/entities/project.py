import typing
from datetime import datetime

import pytz
from sqlalchemy import and_, Column

from .entity import Entity
from ..utils.uuid import generate_id
from ..corpora_orm import DbProject, DbProjectLink, ProjectStatus


class Project(Entity):
    table = DbProject
    list_entities = [DbProject.created_at, DbProject.name, DbProject.id]

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
    ) -> "Project":
        """
        Create a new Project and related objects and store in the database. UUIDs are generated for all new table
        entries.

        """
        uuid = generate_id()

        # Setting Defaults
        links = links if links else []

        #  Prevent accidentally linking an existing row to a different Project. This maintains the relationship of one
        #  to many for links
        [link.pop("id", None) for link in links]  # sanitize of ids

        new_db_object = DbProject(
            id=uuid,
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
                links, DbProjectLink, add_columns=dict(project_id=uuid, project_status=status)
            ),
        )

        cls.db.session.add(new_db_object)
        cls.db.session.commit()
        return cls(new_db_object)

    @classmethod
    def get_project(cls, project_uuid):
        """
        Given the project_uuid, retrieve a live project.
        :param project_uuid:
        :return:
        """
        return cls.get((project_uuid, ProjectStatus.LIVE.name))

    @classmethod
    def list_in_time_range(
        cls, to_date: float = None, from_date: float = None, list_entities: typing.List[Column] = None
    ) -> typing.List[typing.Dict]:
        """

        :param to_date: Filter dates earlier than this. Unix timestamp since the epoch in UTC timezone.
        :param from_date: Filter dates later than this. Unix timestamp since the epoch in UTC timezone.
        :param list_entities: The columns to retrieve from the table.
        :return: The results is a list of flattened dictionaries containing the `list_entities`
        """

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _result[_field] = getattr(db_object, _field)
            return _result

        filters = [DbProject.status == ProjectStatus.LIVE.name]
        list_entities = list_entities if list_entities else cls.list_entities
        if to_date:
            filters.append(DbProject.created_at <= datetime.fromtimestamp(to_date, tz=pytz.UTC))
        if from_date:
            filters.append(DbProject.created_at >= datetime.fromtimestamp(from_date, tz=pytz.UTC))

        results = [
            to_dict(result)
            for result in cls.db.session.query(DbProject).with_entities(*list_entities).filter(and_(*filters)).all()
        ]

        return results
