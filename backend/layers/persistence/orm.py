from sqlalchemy import Column, DateTime, Enum, ForeignKey, String
from sqlalchemy.dialects.postgresql import ARRAY, BOOLEAN, JSON, UUID
from sqlalchemy.orm import registry
from sqlalchemy.schema import MetaData

from backend.layers.common.entities import DatasetArtifactType
from backend.layers.persistence.constants import SCHEMA_NAME

metadata = MetaData(schema=SCHEMA_NAME)
mapper_registry = registry(metadata=metadata)


@mapper_registry.mapped
class CollectionTable:

    __tablename__ = "Collection"

    id = Column(UUID(as_uuid=True), primary_key=True)
    version_id = Column(UUID(as_uuid=True))
    originally_published_at = Column(DateTime)
    revised_at = Column(DateTime)
    tombstone = Column(BOOLEAN)


@mapper_registry.mapped
class CollectionVersionTable:

    __tablename__ = "CollectionVersion"

    id = Column(UUID(as_uuid=True), primary_key=True)
    collection_id = Column(UUID(as_uuid=True))
    collection_metadata = Column(JSON)
    owner = Column(String)
    curator_name = Column(String)
    publisher_metadata = Column(JSON)
    published_at = Column(DateTime)
    created_at = Column(DateTime)
    schema_version = Column(String)
    datasets = Column(ARRAY(UUID(as_uuid=True)))


@mapper_registry.mapped
class DatasetTable:

    __tablename__ = "Dataset"

    id = Column(UUID(as_uuid=True), primary_key=True)
    version_id = Column(UUID(as_uuid=True))
    published_at = Column(DateTime)
    tombstone = Column(BOOLEAN)


@mapper_registry.mapped
class DatasetVersionTable:

    __tablename__ = "DatasetVersion"

    id = Column(UUID(as_uuid=True), primary_key=True)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("Dataset.id"))
    collection_id = Column(UUID(as_uuid=True))
    created_at = Column(DateTime)
    dataset_metadata = Column(JSON)
    artifacts = Column(ARRAY(UUID(as_uuid=True)))
    status = Column(JSON)


@mapper_registry.mapped
class DatasetArtifactTable:

    __tablename__ = "DatasetArtifact"

    id = Column(UUID(as_uuid=True), primary_key=True)
    type = Column(Enum(DatasetArtifactType))
    uri = Column(String)
