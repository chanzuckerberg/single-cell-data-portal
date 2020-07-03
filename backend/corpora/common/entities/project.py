from .entity import Entity
from ..corpora_orm import DbProject


class Project(Entity):
    table = DbProject

    def __init__(self, db_object: DbProject):
        super().__init__(db_object)
