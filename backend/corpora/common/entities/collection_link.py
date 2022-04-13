import typing

from sqlalchemy.orm import Session

from backend.corpora.common.corpora_orm import DbCollectionLink
from backend.corpora.common.entities.entity import Entity


class CollectionLink(Entity):
    table = DbCollectionLink

    @classmethod
    def create(cls, session: Session, collection_id: str, **kwargs) -> "CollectionLink":
        link = DbCollectionLink(collection_id=collection_id, **kwargs)
        session.add(link)
        session.commit()
        return cls(link)

    @classmethod
    def retrieve_all_links_for_a_collection(cls, session: Session, collection_id: str) -> typing.List["CollectionLink"]:
        return session.query(cls.table).filter(cls.table.collection_id == collection_id).all()
