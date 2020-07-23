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
        result["owner"] = result.pop("user")
        result["links"] = [dict(url=link["link_url"], type=link["link_type"]) for link in result["links"]]
        result["attestation"] = dict(needed=result.pop("needs_attestation", None), tc_uri=result.pop("tc_uri", None))
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = dataset.pop("deployment_directories")
            dataset["dataset_assets"] = dataset.pop("artifacts")
            dataset["preprint_doi"] = dict(title=dataset.pop("preprint_doi"))
            dataset["publication_doi"] = dict(title=dataset.pop("publication_doi"))
        return result
