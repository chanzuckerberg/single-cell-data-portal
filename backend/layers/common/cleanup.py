from backend.layers.business.entities import CollectionMetadataUpdate
from backend.layers.common.entities import CollectionMetadata
from typing import Union


def strip_fields(metadata: Union[CollectionMetadata, CollectionMetadataUpdate]):
    if metadata.name is not None:
        metadata.name = metadata.name.strip()
    if metadata.description is not None:
        metadata.description = metadata.description.strip()
    if metadata.contact_name is not None:
        metadata.contact_name = metadata.contact_name.strip()
    if metadata.contact_email is not None:
        metadata.contact_email = metadata.contact_email.strip()
    if metadata.links is not None:
        for link in metadata.links:
            link.strip_fields()
    if metadata.consortia is not None:
        metadata.consortia = [consortium.strip() for consortium in metadata.consortia]


def sort_consortia(metadata: Union[CollectionMetadata, CollectionMetadataUpdate]):
    if metadata.consortia is not None:
        metadata.consortia.sort()


def sanitize(metadata: Union[CollectionMetadata, CollectionMetadataUpdate]):
    strip_fields(metadata)
    sort_consortia(metadata)
