#!/usr/bin/env bash

# Mirror RDS (postgres db) data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging). Steps:
# 1. Dump the database using pg_dump from the src env
# 2. Loads the pg_dump data file into a PostgreSQL db in the dest env
# 3. Updates table columns that contain deployment environment-specific URLs to point at the dest environment
#
# THIS IS DESTRUCTIVE! (It will not run in prod)

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh

echo Mirroring RDS data from $SRC_ENV to $DEST_ENV

kill_ssh_tunnel()
{
  pkill -f "bastion\..*\.single-cell\.czi\.technology" || true
}

DB_DUMP_FILE=`mktemp`

export DEPLOYMENT_STAGE=$SRC_ENV
if [[ $SRC_ENV == 'staging' ]]; then
   export AWS_PROFILE=single-cell-dev
else
   export AWS_PROFILE=single-cell-${SRC_ENV}
fi
# TODO: move kill ssh tunnel to Makefile db/tunnel command
kill_ssh_tunnel
cd $SCRIPTS_DIR/..
make db/tunnel
pgrep -fl bastion
make db/dump OUTFILE=$DB_DUMP_FILE

export DEPLOYMENT_STAGE=$DEST_ENV
export AWS_PROFILE=single-cell-dev
kill_ssh_tunnel
make db/tunnel
pgrep -fl bastion

#  For safety, also dump the destination db to a local file, just in case.
DEST_DB_BACKUP_DUMP_FILE="${DEST_ENV}_"`date +%Y%m%d_%H%M%S`".sqlc"
make db/dump OUTFILE=$DEST_DB_BACKUP_DUMP_FILE
echo Created backup dump of destination database: $DEST_DB_BACKUP_DUMP_FILE

DB_PW=`aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'`

DB_NAME="corpora_${DEPLOYMENT_STAGE}"
DB_USER=corpora_${DEPLOYMENT_STAGE}

PGPASSWORD=${DB_PW} pg_restore --clean --if-exists --no-owner --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} ${DB_DUMP_FILE}

PGPASSWORD=${DB_PW} psql --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} -c "UPDATE dataset SET explorer_url = regexp_replace(explorer_url, '(https:\\/\\/)(.+?)(\\/.+)', '\\1cellxgene.${DEPLOYMENT_STAGE}.single-cell.czi.technology\\3') WHERE explorer_url IS NOT NULL"

PGPASSWORD=${DB_PW} psql --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} -c "UPDATE dataset_artifact SET s3_uri = regexp_replace(s3_uri, '(s3:\\/\\/)([[:alpha:]]+-[[:alpha:]]+-)([[:alpha:]]+)(\\/.+)', '\\1\\2${DEPLOYMENT_STAGE}\\4') WHERE s3_uri IS NOT NULL"

kill_ssh_tunnel

pgrep -fl bastion
