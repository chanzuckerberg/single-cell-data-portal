import csv
import os
import typing
from collections import OrderedDict
from datetime import datetime
from pathlib import PurePosixPath

from ddtrace import tracer
from sqlalchemy import or_
from sqlalchemy.orm import Session

from urllib.parse import urlparse

from backend.common.entities.dataset_asset import DatasetAsset
from backend.common.entities.entity import Entity
from backend.common.entities.geneset import Geneset
from backend.common.corpora_orm import (
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    UploadStatus,
    ProcessingStatus,
    DbGenesetDatasetLink,
    generate_id,
    DatasetArtifactFileType,
    ConversionStatus,
)
from backend.common.utils.corpora_constants import CorporaConstants
from backend.common.utils.db_helpers import clone
from backend.common.utils.ontology_mappings.ontology_map_loader import ontology_mappings
from backend.common.utils.s3_buckets import buckets


class Dataset(Entity):
    table = DbDataset

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)

    @classmethod
    @tracer.wrap()
    def create(
        cls,
        session: Session,
        revision: int = 0,
        name: str = "",
        artifacts: list = None,
        processing_status: dict = None,
        **kwargs,
    ) -> "Dataset":
        """
        Creates a new dataset and related objects and store in the database. IDs are generated for all new table
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
    @tracer.wrap()
    def get(
        cls, session: Session, dataset_id=None, include_tombstones=False, collection_id=None
    ) -> typing.Optional["Dataset"]:
        filters = []
        if not include_tombstones:
            filters.append(cls.table.tombstone != True)  # noqa
        if collection_id:
            filters.append(cls.table.collection_id == collection_id)
        if dataset_id:
            filters.append(cls.table.id == dataset_id)
        result = session.query(cls.table).filter(*filters).one_or_none()
        dataset = cls(result) if result else None
        return dataset

    @classmethod
    @tracer.wrap()
    def get_by_explorer_url(cls, session: Session, explorer_url):
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

    @classmethod
    @tracer.wrap()
    def list(
        cls, session: Session, collection_id: str = None, include_tombstones: bool = False
    ) -> typing.List["Dataset"]:
        """
        Retrieves a list of Datasets from the database

        :param session: an open database session
        :param include_tombstones: True of to include tombstoned datasets
        :return: A list of Datasets
        """
        filters = []

        if not include_tombstones:
            filters.append(or_(cls.table.tombstone == False, cls.table.tombstone == None))  # noqa
        if collection_id:
            filters.append(cls.table.collection_id == collection_id)

        if filters:
            return [cls(obj) for obj in session.query(cls.table).filter(*filters).all()]
        else:
            return [cls(obj) for obj in session.query(cls.table)]

    def get_asset(self, asset_id) -> typing.Union[DatasetAsset, None]:
        """
        Retrieve the asset if it exists in the dataset.
        :param asset_id: uuid of the asset to find
        :return: If the asset is found it is returned, else None is returned.
        """
        asset = [asset for asset in self.artifacts if asset.id == asset_id]
        return None if not asset else DatasetAsset(asset[0])

    def get_assets(self):
        """
        Retrieve all the assets for the dataset
        """
        assets = [
            asset for asset in self.artifacts if asset.filename != CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME
        ]
        return assets

    @staticmethod
    def transform_sex_for_schema_2_0_0(dataset):
        schema_version = dataset.get("schema_version")
        if "sex" in dataset and (schema_version is None or schema_version < "2.0.0"):
            dataset["sex"] = [{"label": s, "sex_ontology_term_id": "unknown"} for s in dataset["sex"]]

    @staticmethod
    def transform_organism_for_schema_2_0_0(dataset):
        schema_version = dataset.get("schema_version")
        # If organism is an object (version 1.1.0), wrap it into an array to be 2.0.0 compliant
        if "organism" in dataset and (schema_version is None or schema_version < "2.0.0"):
            dataset["organism"] = [dataset["organism"]]

    @staticmethod
    def enrich_development_stage_with_ancestors(dataset):
        Dataset._enrich_with_ancestors(
            dataset, "development_stage", ontology_mappings.development_stage_ontology_mapping
        )

    @staticmethod
    def enrich_tissue_with_ancestors(dataset):
        """
        Tag dataset with ancestors for all tissues in the given dataset, if any.
        """
        Dataset._enrich_with_ancestors(dataset, "tissue", ontology_mappings.tissue_ontology_mapping)

    @staticmethod
    def enrich_cell_type_with_ancestors(dataset):
        """
        Tag dataset with ancestors for all cell types in the given dataset, if any.
        """
        Dataset._enrich_with_ancestors(dataset, "cell_type", ontology_mappings.cell_type_ontology_mapping)

    def _enrich_with_ancestors(dataset, key, ontology_mapping):
        """
        Tag dataset with ancestors for all values of the given key, if any.
        """
        if key not in dataset:
            return

        terms = [e["ontology_term_id"] for e in dataset[key]]

        if not terms:
            return

        ancestors = [ontology_mapping.get(term) for term in terms]
        flattened_ancestors = [item for sublist in ancestors if sublist for item in sublist]
        unique_ancestors = list(OrderedDict.fromkeys(flattened_ancestors))
        if unique_ancestors:
            dataset[f"{key}_ancestors"] = unique_ancestors

    def _create_new_explorer_url(self, new_id: str) -> str:
        if self.explorer_url is None:
            return None
        original_url = urlparse(self.explorer_url)
        original_path = PurePosixPath(original_url.path)
        new_path = str(original_path.parent / f"{new_id}.cxg/")
        new_url = original_url._replace(path=new_path).geturl()
        # Note: the final slash is mandatory, otherwise the explorer won't load this link
        return f"{new_url}/"

    def create_revision(self, revision_collection_id: str) -> "Dataset":
        """
        Generate a dataset revision from a dataset in a public collection
        :param revision_collection_id: specify the collection revision to which this dataset revision belongs
        :return: dataset revision.

        """
        revision_dataset_id = generate_id()
        revision_explorer_url = self._create_new_explorer_url(revision_dataset_id)
        revision_dataset = clone(
            self.db_object,
            id=revision_dataset_id,
            collection_id=revision_collection_id,
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
        buckets.explorer_bucket.upload_file(csv_file, s3_file)
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
            self_reported_ethnicity=None,
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
        collection_id = self.collection.revision_of if self.collection.revision_of else self.collection.id
        self.update(collection_id=collection_id, published=True, published_at=now, commit=False)

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
                    if (
                        artifact.filetype == DatasetArtifactFileType.RDS
                        and revision.processing_status.rds_status == ConversionStatus.SKIPPED
                    ):
                        # Delete old .rds (Seurat) dataset artifact clones if rds_status for dataset revision is SKIPPED
                        DatasetAsset(artifact).delete(commit=False)
                    else:
                        artifact.dataset_id = self.id

            elif revision.tombstone:
                # tombstone
                revision.tombstone_dataset_and_delete_child_objects()

            updates = revision.to_dict(
                remove_attr=[
                    "updated_at",
                    "created_at",
                    "id",
                    "collection_id",
                    "original_id",
                    "published",
                    "revised_at",
                    "explorer_url",
                ],
                remove_relationships=True,
            )

            if revision.tombstone is not False:
                self.update(commit=False, **updates)
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
