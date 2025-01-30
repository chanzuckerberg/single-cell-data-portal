# Load RDS (postgres db) data from a (prod) db dump on S3
# to a specified destination deployment environment (dev or staging)
#
# THIS IS DESTRUCTIVE for the destination env! The source env will
# never be modified, but the dest env's data will be replaced.
#
# This script is intended to be run in as an AWS ECS task, where the IAM role is appropriately configured to allow for
# dest env S3 bucket RDS access (dev or staging).

LOCAL_DB_DUMP_FILE=db_dump.psql
if [[ -z $DEPLOYMENT_STAGE ]]; then
    echo "DEPLOYMENT_STAGE is not set"
    exit 1
fi
if [[ $DEPLOYMENT_STAGE -eq "prod" ]]; then
    echo "DEPLOYMENT_STAGE is set to prod, which is not allowed"
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

eval `aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString' | jq -r '.database_uri | match("postgresql://(.+):(.+)@(.+)/(.+)").captures | "DB_USER=\(.[0].string) PGPASSWORD=\(.[1].string) DB_HOST=\(.[2].string) DB_NAME=\(.[3].string)"'`
aws s3 cp $DB_DUMP_S3_URI $LOCAL_DB_DUMP_FILE
# TODO: Remove echo statements once we're sure this is working
echo pg_restore --host $DB_HOST --dbname=$DB_NAME --username $DB_USER --clean --if-exists --no-owner --no-privileges --no-comments --schema=persistence_schema db_dump.psql
echo psql --host $DB_HOST --dbname=$DB_NAME --username $DB_USER -c "UPDATE persistence_schema.\"DatasetArtifact\" SET uri = regexp_replace(uri, '(s3:\\/\\/)([[:alpha:]]+-[[:alpha:]]+-)([[:alpha:]]+)(\\/.+)', '\\1\\2staging\\4') WHERE uri IS NOT NULL;"
