import typing
import uuid

from .dataset_asset import DatasetAsset
from .entity import Entity
from ..corpora_orm import DbDataset, DbDatasetArtifact, DbDeploymentDirectory


class Dataset(Entity):
    table = DbDataset

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        revision: int = 0,
        name: str = "",
        organism: dict = None,
        tissue: list = None,
        assay: list = None,
        disease: list = None,
        sex: list = None,
        ethnicity: list = None,
        development_stage: list = None,
        artifacts: list = None,
        deployment_directories: list = None,
        **kwargs,
    ) -> "Dataset":
        """
        Creates a new dataset and related objects and store in the database. UUIDs are generated for all new table
        entries.
        """
        primary_key = str(uuid.uuid4())

        # Setting Defaults
        artifacts = artifacts if artifacts else []
        deployment_directories = deployment_directories if deployment_directories else []

        new_db_object = DbDataset(
            id=primary_key,
            revision=revision,
            name=name,
            organism=organism,
            tissue=tissue,
            assay=assay,
            disease=disease,
            sex=sex,
            ethnicity=ethnicity,
            development_stage=development_stage,
            artifacts=cls._create_sub_objects(artifacts, DbDatasetArtifact, add_columns=dict(dataset_id=primary_key)),
            deployment_directories=cls._create_sub_objects(
                deployment_directories,
                DbDeploymentDirectory,
                add_columns=dict(dataset_id=primary_key),
            ),
            **kwargs,
        )

        cls.db.session.add(new_db_object)
        cls.db.session.flush()
        cls.db.commit()

        return cls(new_db_object)

    def get_asset(self, asset_uuid) -> typing.Union[DatasetAsset, None]:
        """
        Retrieve the asset if it exists in the dataset.
        :param asset_uuid: uuid of the asset to find
        :return: If the asset is found it is returned, else None is returned.
        """
        asset = [asset for asset in self.artifacts if asset.id == asset_uuid]
        return None if not asset else DatasetAsset(asset[0])

    def delete(self):
        """
        Delete the Dataset and all child objects.
        """
        super().delete()
