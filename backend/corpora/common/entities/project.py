from .entity import Entity
from ..corpora_orm import DbProject, DbProjectLink


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
    ) -> "Project":
        """
        Need to check if one exists before creating
        :param id:
        :param status:
        :param name:
        :param description:
        :param owner:
        :param s3_bucket:
        :param tc_uri:
        :param needs_attestation:
        :param processing_state:
        :param validation_state:
        :param links:
        :return:
        """
        uuid = cls.db.generate_id(DbProject, status)

        """
        Prevent accidentally linking an existing row to a different Project. This maintains the relationship of one
        to many for links
        """
        links = links if links else []
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
            links=cls._create_sub_objects(links, DbProjectLink, project_id=uuid, project_status=status),
        )

        cls.db.session.add(new_db_object)
        cls.db.session.commit()
        return cls(new_db_object)
