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
S3_SYNC_CMD="aws s3 sync --no-progress"

set -x
$S3_SYNC_CMD --exclude '*loom' s3://corpora-data-${SRC_ENV}/ s3://corpora-data-${DEST_ENV}/
$S3_SYNC_CMD s3://hosted-cellxgene-${SRC_ENV}/ s3://hosted-cellxgene-${DEST_ENV}/
$S3_SYNC_CMD s3://cellxgene-wmg-${SRC_ENV}/ s3://cellxgene-wmg-${DEST_ENV}/
set +x

