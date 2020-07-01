import typing

from .entity import Entity
from ..corpora_orm import DbDataset


class Dataset(Entity):
    tables = [
        DbDataset,
    ]

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)

    @classmethod
    def _query(cls, key: str) -> typing.List["Dataset"]:
        """
        Given a key representing a dataset UUID, queries the database to return the associated dataset data.
        :param key: dataset.id is the primary key.
        :return: list of query result rows
        """
        result = cls.db.query(table_args=[DbDataset], filter_args=[key == DbDataset.id])
        return result

    @classmethod
    def _load(cls, db_result: typing.List["Dataset"]) -> "Dataset":
        """
        Parses a database query result into a Dataset Entity instance.
        :param db_result: list of query result rows
        :return: Entity
        """
        try:
            return cls(db_result[0])
        except IndexError:
            return None

    def __getattr__(self, name):
        return self.db_obj.__getattribute__(name)
