import typing

from .dataset_asset import DatasetAsset
from .entity import Entity
from ..corpora_orm import (
    DbDataset,
    DbDatasetArtifact,
    DbDeploymentDirectory,
    DbDatasetProcessingStatus,
    UploadStatus,
    ProcessingStatus,
    DbGenesetDatasetLink,
)


class Dataset(Entity):
    table = DbDataset

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        session,
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
        processing_status: dict = None,
        **kwargs,
    ) -> "Dataset":
        """
        Creates a new dataset and related objects and store in the database. UUIDs are generated for all new table
        entries.
        """
        dataset = DbDataset(
            revision=revision,
            name=name,
            organism=organism,
            tissue=tissue,
            assay=assay,
            disease=disease,
            sex=sex,
            ethnicity=ethnicity,
            development_stage=development_stage,
            **kwargs,
        )
        if artifacts:
            dataset.artifacts = [DbDatasetArtifact(dataset_id=dataset.id, **art) for art in artifacts]
        if deployment_directories:
            dataset.deployment_directories = [
                DbDeploymentDirectory(dataset_id=dataset.id, **dd) for dd in deployment_directories
            ]
        processing_status = processing_status if processing_status else {}
        dataset.processing_status = DbDatasetProcessingStatus(dataset_id=dataset.id, **processing_status)
        session.add(dataset)
        session.commit()

        return cls(dataset)

    def update(
        self, artifacts: list = None, deployment_directories: list = None, processing_status: dict = None, **kwargs
    ) -> None:
        """
        Update an existing dataset to match provided the parameters. The specified column are replaced.
        :param artifacts: Artifacts to create and connect to the dataset. If present, the existing attached entries will
         be removed and replaced with new entries.
        :param deployment_directories: Deployment directories to create and connect to the dataset. If present, the
         existing attached entries will be removed and replaced with new entries.
        :param processing_status: A Processing status entity to create and connect to the dataset. If present, the
         existing attached entries will be removed and replaced with new entries.
        :param kwargs: Any other fields in the dataset that will be replaced.
        """
        if artifacts or deployment_directories or processing_status:
            if artifacts:
                for af in self.artifacts:
                    self.session.delete(af)
                new_objs = [DbDatasetArtifact(dataset_id=self.id, **art) for art in artifacts]
                self.session.add_all(new_objs)
            if deployment_directories:
                for dd in self.deployment_directories:
                    self.session.delete(dd)
                new_objs = [DbDeploymentDirectory(dataset_id=self.id, **dd) for dd in deployment_directories]
                self.session.add_all(new_objs)
            if processing_status:
                if self.processing_status:
                    self.session.delete(self.processing_status)
                new_obj = DbDatasetProcessingStatus(dataset_id=self.id, **processing_status)
                self.session.add(new_obj)

            self.session.flush()

        super().update(**kwargs)
        self.session.commit()

    @classmethod
    def get(cls, session, dataset_uuid, include_tombstones=False):
        dataset = super().get(session, dataset_uuid)
        if not include_tombstones:
            if dataset and dataset.tombstone is True:
                return None
        return dataset

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

    def tombstone_dataset_and_delete_child_objects(self):
        self.update(tombstone=True)
        self.session.delete(self.processing_status)
        for dd in self.deployment_directories:
            self.session.delete(dd)
        for af in self.artifacts:
            self.session.delete(af)
        self.session.query(DbGenesetDatasetLink).filter(DbGenesetDatasetLink.dataset_id == self.id).delete(
            synchronize_session="evaluate"
        )
        self.session.commit()

    def dataset_and_asset_deletion(self):
        for artifact in self.artifacts:
            asset = DatasetAsset.get(self.session, artifact.id)
            asset.delete_from_s3()
        self.tombstone_dataset_and_delete_child_objects()

    @staticmethod
    def new_processing_status():
        return {
            "upload_status": UploadStatus.WAITING,
            "upload_progress": 0,
            "processing_status": ProcessingStatus.PENDING,
        }
