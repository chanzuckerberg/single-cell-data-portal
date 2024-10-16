# Dump RDS (postgres db) data from production env to a db dump file stored on S3.
#
# This script is intended to be run in as an AWS ECS task, where the IAM role is appropriately configured to allow for
# prod env RDS access.

LOCAL_DB_DUMP_FILE=db_dump.psql
if [[ -z $DEPLOYMENT_STAGE ]]; then
    echo "DEPLOYMENT_STAGE is not set"
    exit 1
fi
if [[ -z $DB_DUMP_S3_URI ]]; then
    echo "DB_DUMP_S3_URI is not set"
    exit 1
fi
if [[ -z $GHA_COMMIT ]]; then
    echo "GHA_COMMIT is not set"
    exit 1
fi

eval `aws secretsmanager get-secret-value --secret-id corpora/backend/prod/database --region us-west-2 | jq -r '.SecretString' | jq -r '.database_uri | match("postgresql://(.+):(.+)@(.+)/(.+)").captures | "DB_USER=\(.[0].string) PGPASSWORD=\(.[1].string) DB_HOST=\(.[2].string) DB_NAME=\(.[3].string)"'`
pg_dump -Fc --host $DB_HOST --dbname=$DB_NAME --username $DB_USER --file=$LOCAL_DB_DUMP_FILE
aws s3 cp $LOCAL_DB_DUMP_FILE ${DB_DUMP_S3_URI}.${GHA_COMMIT}
