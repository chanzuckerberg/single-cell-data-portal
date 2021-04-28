import typing
from sqlalchemy.exc import SQLAlchemyError

from ..corpora_orm import DbGeneset, DbGenesetDatasetLink
from .entity import Entity
from ..utils.exceptions import CorporaException


class Geneset(Entity):
    table = DbGeneset

    @classmethod
    def create(cls, session, name: str, description: str, genes: list, dataset_ids: list = [], **kwargs):
        gene_set = DbGeneset(name=name, description=description, genes=genes, **kwargs)
        session.add(gene_set)
        session.commit()
        if dataset_ids:
            for dataset_id in dataset_ids:
                GenesetDatasetLink.create(session, geneset_id=gene_set.id, dataset_id=dataset_id)
        return cls(gene_set)

    @classmethod
    def retrieve_all_genesets_for_a_collection(cls, session, collection_id):
        genesets = session.query(cls.table).filter(cls.table.collection_id == collection_id).all()
        reshaped_genesets = []
        for geneset in genesets:
            reshaped_genesets.append(
                geneset.to_dict(
                    remove_attr=[
                        "genes",
                        "created_at",
                        "updated_at",
                        "collection",
                        "collection_id",
                        "collection_visibility",
                    ]
                )
            )
        return reshaped_genesets

    def convert_geneset_to_gene_dicts(self):
        max_additional_params = 0
        gene_rows = []
        for gene in self.genes:
            gene_row = {
                "GENE_SET_NAME": self.name,
                "GENE_SET_DESCRIPTION": self.description,
                "GENE_SYMBOL": gene["gene_symbol"],
                "GENE_DESCRIPTION": gene["gene_description"],
            }
            for key, value in gene.get("additional_params",dict()).items():
                gene_row[key.upper()] = value
            gene_rows.append(gene_row)

            addit_params_count = len(gene.get("additional_params", "")) // 2
            if addit_params_count > max_additional_params:
                max_additional_params = addit_params_count
        return gene_rows, max_additional_params


class GenesetDatasetLink(Entity):
    table = DbGenesetDatasetLink

    @classmethod
    def create(cls, session, geneset_id: str, dataset_id: str):
        link = DbGenesetDatasetLink(geneset_id=geneset_id, dataset_id=dataset_id)
        session.add(link)
        session.commit()
        return cls(link)

    @classmethod
    def get(cls, session, geneset_id: str, dataset_id: str):
        link = (
            session.query(cls.table)
            .filter(cls.table.geneset_id == geneset_id, cls.table.dataset_id == dataset_id)
            .one_or_none()
        )
        if link:
            return cls(link)
        else:
            return None

    @classmethod
    def update_links_for_a_dataset(
        cls, session, dataset_id: str, add: typing.Optional[list] = None, remove: typing.Optional[list] = None
    ):
        for gene_set_id in remove:
            link = cls.get(session=session, dataset_id=dataset_id, geneset_id=gene_set_id)
            if link is None:
                session.rollback()
                raise CorporaException()
            session.delete(link)
        try:
            for geneset_id in add:
                session.add(DbGenesetDatasetLink(geneset_id=geneset_id, dataset_id=dataset_id))
            session.commit()
        except SQLAlchemyError:
            session.rollback()
            raise CorporaException
