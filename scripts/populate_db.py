#!/usr/bin/env python
import os
import sys

import click

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists

from backend.common.corpora_config import CorporaDbConfig

env = os.environ.get("DEPLOYMENT_STAGE")

# Importing tests.unit overwrites our deployment stage env var.
# So we're putting it back here.
os.environ["DEPLOYMENT_STAGE"] = env


@click.command()
@click.option("--create-schema/--skip-create-schema", default=False, help="Create schema if it doesn't exist")
@click.option("--recreate-db/--skip-recreate-db", default=True, help="Drop and recreate DB tables.")
@click.option("--populate-data/--skip-populate", default=True, help="Add test data to db.")
@click.option("--drop-db/--skip-drop-db", default=False, help="Drop Database.")
def run_db_stuff(create_schema, recreate_db, populate_data, drop_db):
    # Create schema.
    if create_schema:
        print(CorporaDbConfig().database_uri)
        engine = create_engine(CorporaDbConfig().database_uri)
        if not database_exists(engine.url):
            print("Database does not exist, creating database")
            # create_database(engine.url)
        else:
            print("Database already exists")
            exit(1)

    # Drop database
    if drop_db:
        engine = create_engine(CorporaDbConfig().database_uri)
        if database_exists(engine.url):
            print("Database exists, dropping database")
            # drop_database(engine.url)
        else:
            print("Database does not exists")
            exit(1)


if __name__ == "__main__":
    run_db_stuff()
