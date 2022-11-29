from dataclasses import dataclass, field
from typing import List
import uuid
from sqlalchemy import Column, DateTime, Enum, ForeignKey, MetaData, String, Table
from sqlalchemy.dialects.postgresql import ARRAY, BOOLEAN, JSON, TEXT, UUID
from sqlalchemy.orm import registry

from backend.layers.common.entities import (
    CanonicalCollection,
    CanonicalDataset,
    CollectionId,
    CollectionVersion as CollectionVersionModel,
    CollectionVersionId,
    DatasetArtifact as DatasetArtifactModel,
    DatasetArtifactType,
    DatasetArtifactId,
    DatasetId,
    DatasetVersion as DatasetVersionModel,
    DatasetVersionId,
)

metadata_obj = MetaData(schema="persistence_schema")
mapper_registry = registry(metadata=metadata_obj)


@mapper_registry.mapped
class Collection(CanonicalCollection):

    __table__ = Table(
        "Collection",
        mapper_registry.metadata,
        Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("version_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("originally_published_at", DateTime),
        Column("tombstoned", BOOLEAN)
    )


@mapper_registry.mapped
class CollectionVersion(CollectionVersionModel):

    canonical_collection: CanonicalCollection = field(default=None)

    __table__ = Table(
        "CollectionVersion",
        mapper_registry.metadata,
        Column("version_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("collection_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("metadata", JSON),
        Column("owner", String),
        Column("publisher_metadata", JSON),
        Column("published_at", DateTime),
        Column("datasets", ARRAY(UUID(as_uuid=True)))
    )


@mapper_registry.mapped
class Dataset(CanonicalDataset):

    __table__ = Table(
        "Dataset",
        mapper_registry.metadata,
        Column("dataset_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("dataset_version_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("published_at", DateTime)
    )


@mapper_registry.mapped
class DatasetVersion(DatasetVersionModel):

    artifacts: List[DatasetArtifactId] = field(default=list())
    canonical_dataset: CanonicalDataset = field(default=None)

    __table__ = Table(
        "DatasetVersion",
        mapper_registry.metadata,
        Column("version_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("dataset_id", UUID(as_uuid=True), ForeignKey("Dataset.dataset_id"), default=uuid.uuid4),
        Column("collection_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("metadata", JSON),
        Column("artifacts", ARRAY(UUID(as_uuid=True))),
        Column("status", JSON)
    )


@mapper_registry.mapped
class DatasetArtifact(DatasetArtifactModel):

    __table__ = Table(
        "DatasetArtifact",
        mapper_registry.metadata,
        Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("type", Enum(DatasetArtifactType)),
        Column("uri", String)
    )


