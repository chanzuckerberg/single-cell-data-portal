import uuid

from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink


class Project(Entity):
    table = DbProject

    def __init__(self, session, db_object: DbProject):
        super().__init__(session, db_object)

    @classmethod
    def create(
        cls,
        session,
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
                session, links, DbProjectLink, add_columns=dict(project_id=primary_key, project_status=status, **kwargs)
            ),
        )

        session.add(new_db_object)
        session.commit()
        return cls(session, new_db_object)
