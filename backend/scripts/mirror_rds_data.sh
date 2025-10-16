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

SCRIPTS_DIR=$(dirname "$0")
source "${SCRIPTS_DIR}/_mirror_utils.sh"
source "${SCRIPTS_DIR}/set_src_dest_envs.sh"

log "Mirroring RDS data from $SRC_ENV to $DEST_ENV"

DB_DUMP_FILE=$(mktemp)

# Dump source database
export DEPLOYMENT_STAGE=$SRC_ENV
set_aws_profile "$SRC_ENV"
cd $SCRIPTS_DIR/..
make db/dump OUTFILE=$DB_DUMP_FILE PORT=${SRC_PORT}

# Prepare destination environment
export DEPLOYMENT_STAGE=$DEST_ENV
set_aws_profile "$DEST_ENV"

if [[ $DEST_ENV != 'rdev' ]]; then
   #  For safety, dump the destination db to a local file, just in case. Not necessary if destination is rdev.
   DEST_DB_BACKUP_DUMP_FILE="${DEST_ENV}_"`date +%Y%m%d_%H%M%S`".sqlc"
   echo "Backing up the destination database to $DEST_DB_BACKUP_DUMP_FILE. Just in case!"
   make db/dump OUTFILE=$DEST_DB_BACKUP_DUMP_FILE PORT=${DEST_PORT}
fi

DB_PW=$(aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region ${AWS_REGION} | jq -r '.SecretString | match(":([^:]*)@").captures[0].string')

# Get database credentials using shared utility
read DB_NAME DB_USER <<< $(get_db_credentials "$DEST_ENV" "$STACK")

# Confirmation prompt (only if not already confirmed in set_src_dest_envs.sh)
if [[ -z "$NO_PROMPT" ]]; then
    read -n 1 -p "ATTENTION: Proceed to replace the destination database \"${DB_NAME}\"? (Y/n) " ANS
    echo
    [[ $ANS == 'Y' ]] || exit 1
fi

function load_src_dump_to_dest_db() {
  PGPASSWORD=${DB_PW} pg_restore --clean --if-exists --no-owner --no-privileges --no-comments --dbname=${DB_NAME} \
  --host 0.0.0.0 --port ${DEST_PORT} --username ${DB_USER} --schema=persistence_schema ${DB_DUMP_FILE}
}

load_src_dump_to_dest_db || load_src_dump_to_dest_db  # Hack for rdev, where it fails on first and succeeds on retry

if [[ $DEST_ENV != 'rdev' ]]; then
  DB_UPDATE_CMDS=$(cat <<EOF
  -c "UPDATE persistence_schema.\"DatasetArtifact\" \
      SET uri = regexp_replace(uri, '^(s3://[^/]*?)-${SRC_ENV}(/.+)', '\1-${DEST_ENV}\2') \
      WHERE uri LIKE 's3://%-${SRC_ENV}/%';"
EOF
)
else
  rdev_bucket_prefix="env-rdev"
  DB_UPDATE_CMDS=$(cat <<EOF
-c "UPDATE persistence_schema.\"DatasetArtifact\" SET uri = regexp_replace(uri, '(s3:\\/\\/)([^/]+)(\\/.+)', '\\1${rdev_bucket_prefix}-artifacts/${STACK}\\3') WHERE uri ~ 's3://corpora-data'; UPDATE persistence_schema.\"DatasetArtifact\" SET uri = regexp_replace(uri, '(s3:\\/\\/)([^/]+)(\\/.+)', '\\1${rdev_bucket_prefix}-cellxgene/${STACK}\\3') WHERE uri ~ 's3://hosted-cellxgene';"
EOF
)
fi

make db/connect ARGS="${DB_UPDATE_CMDS}"
