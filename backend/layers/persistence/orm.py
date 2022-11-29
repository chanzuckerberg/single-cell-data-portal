from dataclasses import dataclass, field
from typing import List
import uuid
from sqlalchemy import Column, DateTime, Enum, ForeignKey, String, Table
from sqlalchemy.dialects.postgresql import ARRAY, BOOLEAN, JSON, TEXT, UUID
from sqlalchemy.orm import registry
from sqlalchemy.schema import MetaData

from backend.layers.common.entities import DatasetArtifactType

metadata = MetaData(schema="persistence_schema")
mapper_registry = registry(metadata=metadata)

@mapper_registry.mapped
class Collection:

    __table__ = Table(
        "Collection",
        mapper_registry.metadata,
        Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("version_id", UUID(as_uuid=True)),
        Column("originally_published_at", DateTime),
        Column("tombstoned", BOOLEAN)
    )


@mapper_registry.mapped
class CollectionVersion:

    __table__ = Table(
        "CollectionVersion",
        mapper_registry.metadata,
        Column("version_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("collection_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("metadata", JSON),
        Column("owner", String),
        Column("publisher_metadata", JSON),
        Column("published_at", DateTime),
        Column("created_at", DateTime),
        Column("datasets", ARRAY(UUID(as_uuid=True)))
    )


@mapper_registry.mapped
class Dataset:

    __table__ = Table(
        "Dataset",
        mapper_registry.metadata,
        Column("dataset_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("dataset_version_id", UUID(as_uuid=True), default=uuid.uuid4),
        Column("published_at", DateTime)
    )


@mapper_registry.mapped
class DatasetVersion:

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
class DatasetArtifact:

    __table__ = Table(
        "DatasetArtifact",
        mapper_registry.metadata,
        Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        Column("type", Enum(DatasetArtifactType)),
        Column("uri", String)
    )


