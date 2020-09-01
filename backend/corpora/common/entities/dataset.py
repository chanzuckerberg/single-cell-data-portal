import typing
import uuid

from backend.corpora.common.entities.asset import Asset
from .entity import Entity
from ..corpora_orm import DbDataset, DbDatasetArtifact, DbDeploymentDirectory, DbContributor, DbDatasetContributor


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
        source_data_location: str = "",
        preprint_doi: str = "",
        publication_doi: str = "",
        artifacts: list = None,
        contributors: list = None,
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
        contributors = contributors if contributors else []

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
            source_data_location=source_data_location,
            preprint_doi=preprint_doi,
            publication_doi=publication_doi,
            artifacts=cls._create_sub_objects(artifacts, DbDatasetArtifact, add_columns=dict(dataset_id=primary_key)),
            deployment_directories=cls._create_sub_objects(
                deployment_directories,
                DbDeploymentDirectory,
                add_columns=dict(dataset_id=primary_key),
            ),
            **kwargs,
        )

        #  Linking many contributors to many datasets
        contributors = cls._create_sub_objects(contributors, DbContributor)
        contributor_dataset_ids = [
            dict(contributor_id=contributor.id, dataset_id=primary_key) for contributor in contributors
        ]
        dataset_contributor = cls._create_sub_objects(contributor_dataset_ids, DbDatasetContributor)

        cls.db.session.add(new_db_object)
        cls.db.session.add_all(contributors)
        cls.db.session.flush()
        cls.db.session.add_all(dataset_contributor)
        cls.db.commit()

        return cls(new_db_object)

    def get_asset(self, asset_uuid) -> typing.Union[Asset, None]:
        """
        Retrieve the asset if it exists in the dataset.
        :param asset_uuid: uuid of the asset to find
        :return: If the asset is found it is returned, else None is returned.
        """
        asset = [asset for asset in self.artifacts if asset.id == asset_uuid]
        if not asset:
            return
        else:
            return Asset(asset[0])

    def delete(self):
        """
        Delete the Dataset and all child objects. Contributors connect to a dataset are deleted if they are not longer
        connected to any datasets.

        TODO: Must also delete s3 objects related to the asset.
        :return:
        """
        contributors = self.db_object.contributors
        for contributor in contributors:
            if len(contributor.datasets) == 1:
                self.db.delete(contributor)
        super().delete()
