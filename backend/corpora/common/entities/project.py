import uuid

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
