#!/usr/bin/env bash

# Sync data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)

set -e

SRC_ENV=${1:-prod}
DEST_ENV=${2:-dev}

if [[ ! "$SRC_ENV" =~ ^(dev|staging|prod)$ ]]; then
  echo ERROR: invalid value for source deployment env: \"$SRC_ENV\"
  exit 1
fi

if [[ ! "$DEST_ENV" =~ ^(dev|staging)$ ]]; then
  echo ERROR: invalid value for destination deployment env: \"$DEST_ENV\"
  exit 1
fi

kill_ssh_tunnel()
{
  pkill -f "bastion\..*\.single-cell\.czi\.technology" || true
}     

# TODO: Remove `--copy-props metadata-directive` once IAM roles have {Get,Put}ObjectTagging perms added
# TODO: Remove --dryrun after testing
S3_SYNC_CMD="aws s3 sync --delete --metadata-directive REPLACE --no-progress --dryrun"

# TODO: Uncomment after testing with below commands
# $S3_SYNC_CMD s3://corpora-data-${SRC_ENV}/ s3://corpora-data-${DEST_ENV}/
# $S3_SYNC_CMD s3://hosted-cellxgene-${SRC_ENV}/ s3://hosted-cellxgene-${DEST_ENV}/
# $S3_SYNC_CMD s3://cellxgene-wmg-${SRC_ENV}/ s3://cellxgene-wmg-${DEST_ENV}/
# $S3_SYNC_CMD s3://corpora-data-${SRC_ENV}/0b696cf4-513e-4e59-b6ac-78d76409e6f8/ s3://atolopko-tmp/corpora-data-${DEST_ENV}/
# $S3_SYNC_CMD s3://hosted-cellxgene-${SRC_ENV}/00099d5e-154f-4a7a-aa8d-fa30c8c0c43c.cxg/ s3://atolopko-tmp/hosted-cellxgene-${DEST_ENV}/
# $S3_SYNC_CMD s3://cellxgene-wmg-${SRC_ENV}/1651599970/ s3://atolopko-tmp/cellxgene-wmg-${DEST_ENV}/


DB_DUMP_FILE=`mktemp`

cd `dirname $0`/..

export DEPLOYMENT_STAGE=$SRC_ENV
if [[ $SRC_ENV == 'staging' ]]; then
   export AWS_PROFILE=single-cell-dev
else
   export AWS_PROFILE=single-cell-${SRC_ENV}
fi
kill_ssh_tunnel
make db/tunnel
pgrep -fl bastion
make db/dump OUTFILE=$DB_DUMP_FILE

export DEPLOYMENT_STAGE=$DEST_ENV
export AWS_PROFILE=single-cell-dev
kill_ssh_tunnel
make db/tunnel
pgrep -fl bastion
./scripts/db_load.sh $DB_DUMP_FILE

kill_ssh_tunnel

pgrep -fl bastion
