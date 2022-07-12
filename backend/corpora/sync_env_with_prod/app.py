import json
import boto3
import subprocess


def copy_relational_db(event, context):
    print("Starting copy")
    try:
        db_uri = generate_db_uri()
        copy_rds_data(db_uri)
        # todo update to write to prod after confirming code works in dev
        write_to_s3("cellxgene-db-dump-dev")
    except Exception as e:
        print(f"Error copying rds data: {e}")
    return {"statusCode": 200, "body": json.dumps("Hello from Lambda!")}


def write_to_s3(bucket_name):
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.upload_file("/tmp/prod_data.sql", "prod_data.sql")


def copy_rds_data(db_uri):
    db_file = "/tmp/prod_data.sql"
    # todo refactor
    db_password = db_uri.split(":")[2].split("@")[0]
    # todo swap
    # dump_command = f"PGPASSWORD={db_password} pg_dump -Fc --dbname=corpora_prod --file={db_file} --host 0.0.0.0 --username corpora_prod" # noqa E501
    dump_command = f"PGPASSWORD={db_password} pg_dump -Fc --dbname=corpora_dev --file={db_file} --host 0.0.0.0 --username corpora_dev"  # noqa E501

    subprocess.Popen(dump_command, shell=True)


def generate_db_uri():
    client = boto3.client("secretsmanager")
    # todo update env to prod
    response = client.get_secret_value(
        SecretId="corpora/backend/dev/database",
    )
    secret = json.loads(response["SecretString"])
    db_uri = secret["database_uri"]
    return db_uri
