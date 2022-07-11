import json
import boto3
import subprocess


def copy_relational_db(event, context):
    print('hey mads')
    db_uri = generate_db_uri()
    copy_rds_data(db_uri)
    write_to_s3('cellxgene-db-dump-prod')
    return {
        'statusCode': 200,
        'body': json.dumps('Hello from Lambda!')
    }


def write_to_s3(bucket_name):
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    bucket.upload_file('/tmp/prod_data.sql', 'prod_data.sql')


def copy_rds_data(db_uri):
    db_file = "/tmp/prod_data.sql"
    # todo refactor
    db_password = db_uri.split(":")[2].split('@')[0]
    dump_command = f"PGPASSWORD={db_password} pg_dump -Fc --dbname=corpora_prod --file={db_file} --host 0.0.0.0 --username corpora_prod"
    subprocess.Popen(dump_command, shell=True)


def generate_db_uri():
    client = boto3.client('secretsmanager')
    response = client.get_secret_value(
        SecretId='corpora/backend/prod/database',
    )
    secret = json.loads(response['SecretString'])
    db_uri = secret['database_uri']
    return db_uri