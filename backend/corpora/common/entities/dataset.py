from .entity import Entity
from ..corpora_orm import DbDataset


class Dataset(Entity):
    table = DbDataset

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)
