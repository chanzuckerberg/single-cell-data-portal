#!/usr/bin/env bash

# Sync S3 data from production environment to a specified destination deployment environment

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

# TODO: Remove `--copy-props metadata-directive` once IAM roles have {Get,Put}ObjectTagging perms added
# TODO: Remove --dryrun after testing
S3_SYNC_CMD="aws s3 sync --delete --copy-props metadata-directive --no-progress --dryrun"


# # TODO: Uncomment after testing with below commands
# # $S3_SYNC_CMD s3://corpora-data-prod/ s3://corpora-data-${DEST_ENV}/
# # $S3_SYNC_CMD s3://hosted-cellxgene-prod/ s3://hosted-cellxgene-${DEST_ENV}/
# # $S3_SYNC_CMD s3://cellxgene-wmg-prod/ s3://cellxgene-wmg-${DEST_ENV}/
# $S3_SYNC_CMD s3://corpora-data-prod/0b696cf4-513e-4e59-b6ac-78d76409e6f8/ s3://atolopko-tmp/corpora-data-${DEST_ENV}/
# $S3_SYNC_CMD s3://hosted-cellxgene-prod/00099d5e-154f-4a7a-aa8d-fa30c8c0c43c.cxg/ s3://atolopko-tmp/hosted-cellxgene-${DEST_ENV}/
# $S3_SYNC_CMD s3://cellxgene-wmg-prod/1651599970/ s3://atolopko-tmp/cellxgene-wmg-${DEST_ENV}/


DB_DUMP_FILE=`mktemp`

cd `dirname $0`/..

export DEPLOYMENT_STAGE=$SRC_ENV
if [[ $SRC_ENV == 'staging' ]]; then
   export AWS_PROFILE=single-cell-dev
else
   export AWS_PROFILE=single-cell-${SRC_ENV}
fi
make db/tunnel
make db/dump OUTFILE=$DB_DUMP_FILE

export DEPLOYMENT_STAGE=$DEST_ENV
export AWS_PROFILE=single-cell-dev
make db/tunnel
./db_load.sh $DB_DUMP_FILE



