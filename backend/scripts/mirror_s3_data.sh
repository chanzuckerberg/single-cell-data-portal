#!/usr/bin/env bash

# Mirror S3 data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)
#
# This will *add* S3 objects to the dest env, but will not remove existing objects (see NOTE, below).
# The src S3 buckets will never be modified.

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh

echo Mirroring S3 data from $SRC_ENV to $DEST_ENV

# TODO: Add --delete once we confirm that is no data in the folders
# that needs to be kept around (buckets are versioned, so we're not in
# jeopardy of losing anything permanently). This would make the
# operation more destructive, of course!
#
# Note: "--copy-props metadata-directive" copies s3 object metadata, but not tags.
# We have not granted s3:GetObjectTagging perm to dev AWS account on prod buckets,
# so this avoids errors. As none of the s3 objects contain tags, this is acceptable.
S3_SYNC_CMD="aws s3 sync --copy-props metadata-directive --no-progress"

set -x
$S3_SYNC_CMD --exclude '*loom' s3://corpora-data-${SRC_ENV}/ s3://corpora-data-${DEST_ENV}/
$S3_SYNC_CMD s3://hosted-cellxgene-${SRC_ENV}/ s3://hosted-cellxgene-${DEST_ENV}/
# We need to assume the sync-datasets-{dev|staging} role, which has the necessary permissions to perform a s3 sync from
# {dev|staging|prod} to {dev|staging} **for the public-access datasets bucket**.
ACCOUNT_ID=`aws sts get-caller-identity | jq -r '.Arn' | grep -o '[0-9]*'`  # First get AWS account id for dev account
SYNC_ROLE_CREDENTIALS=`aws sts assume-role --role-arn arn:aws:iam::${ACCOUNT_ID}:role/sync-datasets-dev \
                                           --role-session-name sync-datasets \
                       | jq -r '.Credentials'`
# Execute sync for datasets bucket with sync role credentials
AWS_SESSION_TOKEN=`echo $SYNC_ROLE_CREDENTIALS | jq -r '.SessionToken'` \
AWS_SECRET_ACCESS_KEY=`echo $SYNC_ROLE_CREDENTIALS | jq -r '.SecretAccessKey'` \
AWS_ACCESS_KEY_ID=`echo $SYNC_ROLE_CREDENTIALS | jq -r '.AccessKeyId'` \
$S3_SYNC_CMD s3://dataset-assets-public-${SRC_ENV}/ s3://dataset-assets-public-${DEST_ENV}/
set +x
