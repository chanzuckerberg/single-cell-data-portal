"""
Drops and recreates all tables according to corpora_orm.py
"""

import os
import sys

from sqlalchemy import create_engine
#
# pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
# sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.corpora_config import CorporaDbConfig
from backend.corpora.common.corpora_orm import Base


from sqlalchemy.schema import DropTable
from sqlalchemy.ext.compiler import compiles
from server.db.cellxgene_orm import Base as CellxGeneBase



@compiles(DropTable, "postgresql")
def _compile_drop_table(element, compiler, **kwargs):
    return compiler.visit_drop_table(element) + " CASCADE"


def create_db():
    engine = create_engine(CorporaDbConfig().database_uri)
    print("Dropping tables")

    Base.metadata.drop_all(engine)
    CellxGeneBase.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)
    CellxGeneBase.metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
