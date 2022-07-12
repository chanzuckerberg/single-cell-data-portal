from backend.corpora.common.corpora_config import CorporaDbConfig
from backend.corpora.sync_env_with_prod.dump_db import copy_relational_db
from backend.corpora.sync_env_with_prod.load_db import load_db


def run(event, context):
    action = event.get("action")
    if action == "dump":
        copy_relational_db(event, context)
    elif action == "load":
        load_db(CorporaDbConfig().database_uri, event.get("db_dump_s3_bucket"), event.get("db_dump_s3_key"))
    else:
        return {"statusCode": 400}

    return {"statusCode": 200}
