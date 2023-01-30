#!/usr/bin/env bash

# Mirror RDS (postgres db) data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging). Steps:
# 1. Dump the database using pg_dump from the src env
# 2. Loads the pg_dump data file into a PostgreSQL db in the dest env
# 3. Updates table columns that contain deployment environment-specific URLs to point at the dest environment
#
# THIS IS DESTRUCTIVE for the dest env! The src env database will never be modified, but dest environment database will
# be replaced. The destination database will be dumped to a local file, for peace of mind.

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh

echo Mirroring RDS data from $SRC_ENV to $DEST_ENV

DB_DUMP_FILE=`mktemp`

export DEPLOYMENT_STAGE=$SRC_ENV
if [[ $SRC_ENV == 'staging' ]]; then
   export AWS_PROFILE=single-cell-dev
else
   export AWS_PROFILE=single-cell-${SRC_ENV}
fi
cd $SCRIPTS_DIR/..
make db/dump OUTFILE=$DB_DUMP_FILE

export DEPLOYMENT_STAGE=$DEST_ENV
export AWS_PROFILE=single-cell-dev

#  For safety, also dump the destination db to a local file, just in case.
DEST_DB_BACKUP_DUMP_FILE="${DEST_ENV}_"`date +%Y%m%d_%H%M%S`".sqlc"
echo "Backing up the destination database to $DEST_DB_BACKUP_DUMP_FILE. Just in case!"
make db/dump OUTFILE=$DEST_DB_BACKUP_DUMP_FILE

DB_PW=`aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'`

DB_NAME="corpora_${DEPLOYMENT_STAGE}"
DB_USER=corpora_${DEPLOYMENT_STAGE}

read -n 1 -p "ATTENTION: Proceed to replace the destination database \"${DB_NAME}\"? (Y/n) " ANS
echo
[[ $ANS == 'Y' ]] || exit 1

make db/tunnel/up
PGPASSWORD=${DB_PW} pg_restore --clean --if-exists --no-owner --no-privileges --no-comments --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} --schema=persistence_schema ${DB_DUMP_FILE}
make db/tunnel/down

DB_UPDATE_CMDS=$(cat <<EOF
-c "UPDATE persistence_schema.\"DatasetArtifact\" SET uri = regexp_replace(uri, '(s3:\\/\\/)([[:alpha:]]+-[[:alpha:]]+-)([[:alpha:]]+)(\\/.+)', '\\1\\2${DEPLOYMENT_STAGE}\\4') WHERE uri IS NOT NULL;"
EOF
)
                 
make db/connect ARGS="${DB_UPDATE_CMDS}"

