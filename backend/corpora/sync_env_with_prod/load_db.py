import os
import re
import subprocess
import tempfile

import boto3
import psycopg2
from sqlalchemy import create_engine

from backend.corpora.common.corpora_orm import Base


def load_db(db_uri, db_dump_s3_bucket, db_dump_s3_key):
    local_prod_data_sql_file = download_db_dump_from_s3(db_dump_s3_bucket, db_dump_s3_key)

    load_rds_data(db_uri, local_prod_data_sql_file)

    env = os.environ.get("DEPLOYMENT_STAGE", "test")
    update_env_specific_urls(db_uri, explorer_host=f"cellxgene.{env}.single-cell.czi.technology", deploy_env=env)
    print("Successfully loaded database with prod data")


def download_db_dump_from_s3(bucket, key):
    s3 = boto3.client("s3")
    _, local_file = tempfile.mkstemp()

    print(f"Downloading prod db dump file from s3://{bucket}/{key}")
    s3.download_file(bucket, key, local_file)

    return local_file


def load_rds_data(db_uri, prod_sql_file):
    db_cxn_args = parse_db_connection_args(db_uri)
    db_name = db_cxn_args["dbname"]

    engine = create_engine(db_uri)

    print("Dropping all tables in database {db_name}...")
    Base.metadata.drop_all(engine)

    # TODO: test in docker psql
    print("Loading prod db dump data into database {db_name}...")
    load_command = (
        f"PGPASSWORD={db_cxn_args['password']} "
        f"psql "
        f"--dbname={db_name} "
        f"--host {db_cxn_args['host']} "
        f"--username {db_cxn_args['user']} "
        f"--file={prod_sql_file}"
    )
    subprocess.run(load_command, shell=True)

    print(f"Done loading database {db_name} with prod db dump data")


def update_env_specific_urls(db_uri, explorer_host, deploy_env):
    db_cxn_args = parse_db_connection_args(db_uri)
    with psycopg2.connect(**db_cxn_args) as conn:
        with conn.cursor() as cur:
            print("Updating dataset.explorer_url values")
            cur.execute(
                f"""
                UPDATE DATASET
                SET explorer_url = regexp_replace(explorer_url, '(https:\\/\\/)(.+?)(\\/.+)', '\\1{explorer_host}\\3')
                WHERE explorer_url IS NOT NULL
                """
            )

            print("Updating dataset_artifact.s3_uri values")
            # noqa: W291
            cur.execute(
                f"""
                UPDATE dataset_artifact
                SET s3_uri =
                  regexp_replace(s3_uri,
                                 '(s3:\\/\\/)([[:alpha:]]+-[[:alpha:]]+-)([[:alpha:]]+)(\\/.+)',
                                 '\\1\\2{deploy_env}\\4')
                WHERE s3_uri IS NOT NULL
                """
            )


def parse_db_connection_args(db_uri):
    patt = re.compile("postgresql://(?P<user>.+):(?P<password>.+)@(?P<host>.+?)(:(?P<port>.+))?/(?P<dbname>.+)")
    match = patt.match(db_uri)
    return dict(
        user=match.group("user"),
        password=match.group("password"),
        host=match.group("host") or 5432,
        port=match.group("port"),
        dbname=match.group("dbname"),
    )


if __name__ == "__main__":
    db_uri = "postgresql://corpora:test_pw@localhost:5432/corpora"
    load_db(db_uri, "cellxgene-db-dump-prod", "data-portal-rds-dump/rds-prod-dump.sql")
