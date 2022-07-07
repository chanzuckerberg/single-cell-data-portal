#!/usr/bin/env bash

# Sync S3 data from production environment to a specified destination deployment environment

set -e

DEST_DEPLOY_ENV=${1:?"expected arg for destination deployment env"}

if [[ ! "$DEST_DEPLOY_ENV" =~ ^(dev|staging)$ ]]; then
  echo invalid value for destination deployment env: \"$DEST_DEPLOY_ENV\"
  exit 1
fi

# TODO: Remove `--copy-props metadata-directive` once IAM roles have {Get,Put}ObjectTagging perms added
S3_SYNC_CMD="aws s3 sync --copy-props metadata-directive --no-progress"

# TODO: Read full source buckets and use real dest buckets
$S3_SYNC_CMD s3://corpora-data-prod/0b696cf4-513e-4e59-b6ac-78d76409e6f8/ s3://atolopko-tmp/corpora-data-${DEST_DEPLOY_ENV}/
$S3_SYNC_CMD s3://hosted-cellxgene-prod/00099d5e-154f-4a7a-aa8d-fa30c8c0c43c.cxg/ s3://atolopko-tmp/hosted-cellxgene-${DEST_DEPLOY_ENV}/
$S3_SYNC_CMD s3://cellxgene-wmg-prod/1651599970/ s3://atolopko-tmp/cellxgene-wmg-${DEST_DEPLOY_ENV}/
