from backend.corpora.common.corpora_orm import DbGeneset, DBGenesetDatasetLink
from backend.corpora.common.entities.entity import Entity


class Geneset(Entity):
    table = DbGeneset

    @classmethod
    def create(cls, session, name: str, description: str, gene_symbols: list, **kwargs):
        gene_set = DbGeneset(name=name, description=description, gene_symbols=gene_symbols, **kwargs)
        session.add(gene_set)
        session.commit()
        return cls(gene_set)


class GenesetDatasetLink(Entity):
    table = DBGenesetDatasetLink

    @classmethod
    def create(cls, session, geneset_id: str, dataset_id: str):
        link = DBGenesetDatasetLink(geneset_id=geneset_id, dataset_id=dataset_id)
        session.add(link)
        session.commit()
        return cls(link)
