from flask import g

from backend.corpora.common.utils.db_session import DBSessionMaker


def get_db():
    if not g.get("db"):
        g.db = DBSessionMaker().session()
    return g.db
