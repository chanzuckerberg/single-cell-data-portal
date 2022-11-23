from dataclasses import dataclass, field
from sqlalchemy import Column, DateTime, Enum, String, Table
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

mapper_registry = registry()


@mapper_registry.mapped
class Collection(CanonicalCollection):

    __table__ = Table(
        "Collection",
        mapper_registry.metadata,
        Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("version_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("originally_published_at", Column(DateTime)),
        Column("tombstoned", Column(BOOLEAN))
    )


@mapper_registry.mapped
class CollectionVersion(CollectionVersionModel):

    canonical_collection: CanonicalCollection = field(default=None)

    __table__ = Table(
        "CollectionVersion",
        mapper_registry.metadata,
        Column("version_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("collection_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("metadata", Column(JSON)),
        Column("owner", Column(String)),
        Column("publisher_metadata", Column(JSON)),
        Column("published_at", Column(DateTime)),
        Column("datasets", Column(ARRAY(UUID(as_uuid=True))))
    )


@mapper_registry.mapped
class Dataset(CanonicalDataset):

    __table__ = Table(
        "Dataset",
        mapper_registry.metadata,
        Column("dataset_id", Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)),
        Column("dataset_version_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
        Column("published_at", Column(DateTime))
    )


@mapper_registry.mapped
class DatasetVersion(DatasetVersionModel):

    artifacts: List[DatasetArtifactId] = field(default=list())
    canonical_dataset: CanonicalDataset = field(default=None)

    __table__ = Table(
        "DatasetVersion",
        mapper_registry.metadata,
        Column("version_id", Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)),
        Column("dataset_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
        Column("collection_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
        Column("metadata", Column(JSON)),
        Column("artifacts", Column(ARRAY(UUID(as_uuid=True)))),
        Column("status", Column(JSON))
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


