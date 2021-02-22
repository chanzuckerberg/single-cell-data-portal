from backend.corpora.common.corpora_orm import DbGeneset
from backend.corpora.common.entities.entity import Entity


class Geneset(Entity):
    table = DbGeneset

    @classmethod
    def create(cls, session, name: str, description: str, gene_symbols: list, **kwargs):
        gene_set = DbGeneset(name=name, description=description, gene_symbols=gene_symbols, **kwargs)
        session.add(gene_set)
        session.commit()
        return cls(gene_set)
