from urllib.parse import urlparse

import csv
import logging
import os
import typing
from datetime import datetime
from pathlib import PurePosixPath

from .dataset_asset import DatasetAsset
from .entity import Entity
from .geneset import Geneset
from ..corpora_orm import (
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    UploadStatus,
    ProcessingStatus,
    DbGenesetDatasetLink,
    CollectionVisibility,
    generate_uuid,
    DatasetArtifactFileType,
)
from ..utils.db_helpers import clone
from ..utils.s3_buckets import cxg_bucket

logger = logging.getLogger(__name__)


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
        artifacts: list = None,
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
            **kwargs,
        )
        if artifacts:
            dataset.artifacts = [DbDatasetArtifact(dataset_id=dataset.id, **art) for art in artifacts]
        processing_status = processing_status if processing_status else {}
        dataset.processing_status = DbDatasetProcessingStatus(dataset_id=dataset.id, **processing_status)
        session.add(dataset)
        session.commit()

        return cls(dataset)

    def update(self, artifacts: list = None, processing_status: dict = None, commit=True, **kwargs) -> None:
        """
        Update an existing dataset to match provided the parameters. The specified column are replaced.
        :param artifacts: Artifacts to create and connect to the dataset. If present, the existing attached entries will
         be removed and replaced with new entries. If an empty list is provided, all dataset artifacts will be deleted.
        :param processing_status: A Processing status entity to create and connect to the dataset. If present, the
         existing attached entries will be removed and replaced with new entries. If an empty dictionary is provided the
         processing_status will be deleted.
        :param kwargs: Any other fields in the dataset that will be replaced.
        """
        if any([i is not None for i in [artifacts, processing_status]]):
            if artifacts is not None:
                for af in self.artifacts:
                    self.session.delete(af)
                new_objs = [DbDatasetArtifact(dataset_id=self.id, **art) for art in artifacts]
                self.session.add_all(new_objs)
            if processing_status is not None:
                if self.processing_status:
                    self.session.delete(self.processing_status)
                if processing_status:
                    new_obj = DbDatasetProcessingStatus(dataset_id=self.id, **processing_status)
                    self.session.add(new_obj)

        super().update(commit=commit, **kwargs)

    @classmethod
    def get(cls, session, dataset_uuid, include_tombstones=False) -> "Dataset":
        dataset = super().get(session, dataset_uuid)
        if not include_tombstones:
            if dataset and dataset.tombstone is True:
                return None
        return dataset

    @classmethod
    def get_by_explorer_url(cls, session, explorer_url):
        """
        Return the most recently created dataset with the given explorer_url or None
        """

        def _get_by_explorer_url_query(url):
            filters = [cls.table.explorer_url == url]
            dataset = session.query(cls.table).filter(*filters).order_by(cls.table.created_at.desc()).limit(1).all()
            return cls(dataset[0]) if dataset else None

        # If dataset cannot be found, look for the url with (or without) the final slash
        dataset = _get_by_explorer_url_query(explorer_url)
        if not dataset:
            if explorer_url[-1] == "/":
                dataset = _get_by_explorer_url_query(explorer_url[:-1])
            else:
                dataset = _get_by_explorer_url_query(f"{explorer_url}/")

        return dataset

    def get_asset(self, asset_uuid) -> typing.Union[DatasetAsset, None]:
        """
        Retrieve the asset if it exists in the dataset.
        :param asset_uuid: uuid of the asset to find
        :return: If the asset is found it is returned, else None is returned.
        """
        asset = [asset for asset in self.artifacts if asset.id == asset_uuid]
        return None if not asset else DatasetAsset(asset[0])

    def _create_new_explorer_url(self, new_uuid: str) -> str:
        if self.explorer_url is None:
            return None
        original_url = urlparse(self.explorer_url)
        original_path = PurePosixPath(original_url.path)
        new_path = str(original_path.parent / f"{new_uuid}.cxg/")
        return original_url._replace(path=new_path).geturl()

    def create_revision(self) -> "Dataset":
        """
        Generate a dataset revision from a dataset in a public collection
        :param dataset_id: specify the uuid if the revision, otherwise one is generated.
        :return: dataset revision.

        """
        revision_dataset_uuid = generate_uuid()
        revision_explorer_url = self._create_new_explorer_url(revision_dataset_uuid)
        revision_dataset = clone(
            self.db_object,
            id=revision_dataset_uuid,
            collection_id=self.collection_id,
            collection_visibility=CollectionVisibility.PRIVATE,
            original_id=self.id,
            explorer_url=revision_explorer_url,
        )
        self.session.add(revision_dataset)
        for artifact in self.artifacts:
            self.session.add(clone(artifact, dataset_id=revision_dataset.id))
        if self.processing_status:
            self.session.add(clone(self.processing_status, dataset_id=revision_dataset.id))
        self.session.commit()
        return Dataset(revision_dataset)

    def tombstone_dataset_and_delete_child_objects(self):
        self.update(tombstone=True, artifacts=[], processing_status={})
        self.session.query(DbGenesetDatasetLink).filter(DbGenesetDatasetLink.dataset_id == self.id).delete(
            synchronize_session="evaluate"
        )

    def asset_deletion(self):
        for artifact in self.artifacts:
            asset = DatasetAsset.get(self.session, artifact.id)
            asset.queue_s3_asset_for_deletion()
            asset.delete(commit=False)

    @staticmethod
    def new_processing_status() -> dict:
        return {
            "upload_status": UploadStatus.WAITING,
            "upload_progress": 0,
            "processing_status": ProcessingStatus.PENDING,
        }

    def copy_csv_to_s3(self, csv_file: str) -> str:
        object_name = get_cxg_bucket_path(self.explorer_url)
        s3_file = f"{object_name}-genesets.csv"
        cxg_bucket.upload_file(csv_file, s3_file)
        return s3_file

    def generate_tidy_csv_for_all_linked_genesets(self, csv_file_path: str) -> str:
        csv_file = os.path.join(csv_file_path, "geneset.csv")
        fieldnames = ["GENE_SET_NAME", "GENE_SET_DESCRIPTION", "GENE_SYMBOL", "GENE_DESCRIPTION"]
        genesets = []
        max_additional_params = 0
        for geneset in self.genesets:
            geneset_entity = Geneset(geneset)
            gene_rows, gene_max = geneset_entity.convert_geneset_to_gene_dicts()
            if gene_max > max_additional_params:
                max_additional_params = gene_max
            genesets += gene_rows

        for i in range(1, max_additional_params + 1):
            fieldnames.append(f"PROVENANCE{i}")
            fieldnames.append(f"PROVENANCE{i}_DESCRIPTION")

        with open(csv_file, "w") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for gene in genesets:
                writer.writerow(gene)
        return csv_file

    def reprocess(self):
        if not self.published:
            self.asset_deletion()
        self.update(
            name="",
            organism=None,
            tissue=None,
            assay=None,
            disease=None,
            sex=None,
            ethnicity=None,
            development_stage=None,
            cell_type=None,
            published=False,
            revision=self.revision + 1,
            explorer_url=None,
            artifacts=[],
        )

    def get_most_recent_artifact(self, filetype=DatasetArtifactFileType.CXG):
        filters = [DbDatasetArtifact.dataset_id == self.id, DbDatasetArtifact.filetype == filetype]
        artifact = (
            self.session.query(DbDatasetArtifact)
            .filter(*filters)
            .order_by(DbDatasetArtifact.created_at.desc())
            .limit(1)
            .all()
        )
        return DatasetAsset(artifact[0]) if artifact else None

    def publish_new(self, now: datetime):
        """
        Publish a new dataset with the published_at datetime populated
        with the provided datetime.
        :param now: Datetime to populate dataset's published_at.
        """
        self.update(collection_visibility=CollectionVisibility.PUBLIC, published=True, published_at=now, commit=False)

    def publish_revision(self, revision: "Dataset", now: datetime) -> bool:
        """
        Publish a revision of a dataset if the dataset under revision is
        different from the existing dataset.
        :return: True if revision differs from existing dataset, else False.
        """
        if revision.tombstone or revision.revision > self.revision:
            # If the revision is different from the original
            self.asset_deletion()
            if revision.revision > self.revision:
                # connect revised artifacts with published dataset
                for artifact in revision.artifacts:
                    artifact.dataset_id = self.id
            elif revision.tombstone:
                # tombstone
                revision.tombstone_dataset_and_delete_child_objects()

            updates = revision.to_dict(
                remove_attr=[
                    "updated_at",
                    "created_at",
                    "collection_visibility",
                    "id",
                    "original_id",
                    "published",
                    "revised_at",
                    "explorer_url",
                ],
                remove_relationships=True,
            )

            if revision.tombstone is not False:
                self.update(commit=False, **updates, remove_attr="published_at")
            else:
                # There was an update to a dataset, so update revised_at
                self.update(commit=False, **updates, revised_at=now)

            return True

        return False


def get_cxg_bucket_path(explorer_url: str) -> str:
    """Parses the S3 cellxgene bucket object prefix for all resources related to this dataset from the explorer_url"""
    object_name = urlparse(explorer_url).path.split("/", 2)[2]
    if object_name.endswith("/"):
        object_name = object_name[:-1]
    base, _ = os.path.splitext(object_name)
    return base
