import typing

from ddtrace import tracer
from sqlalchemy.orm import Session

from backend.common.corpora_orm import DbCollectionLink
from backend.common.entities.entity import Entity


class CollectionLink(Entity):
    table = DbCollectionLink

    @classmethod
    @tracer.wrap()
    def create(cls, session: Session, collection_id: str, **kwargs) -> "CollectionLink":
        link = DbCollectionLink(collection_id=collection_id, **kwargs)
        session.add(link)
        session.commit()
        return cls(link)

    @classmethod
    @tracer.wrap()
    def retrieve_all_links_for_a_collection(cls, session: Session, collection_id: str) -> typing.List["CollectionLink"]:
        return session.query(cls.table).filter(cls.table.collection_id == collection_id).all()
