#!/usr/bin/env bash

# Mirror S3 data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)
#
# This will *add* S3 objects to the dest env, but will not remove existing objects (see NOTE, below).
# The src S3 buckets will never be modified.
#
# This script is intended to be run in as an AWS ECS task, where the IAM role is appropriately configured to allow for
# prod env S3 bucket reads and dest env S3 bucket writes (dev or staging).

set -e

if [[ -z $DEPLOYMENT_STAGE ]]; then
    echo "DEPLOYMENT_STAGE is not set"
    exit 1
fi

if [[ $DEPLOYMENT_STAGE -eq "prod" ]]; then
    echo "DEPLOYMENT_STAGE is set to prod, which is not allowed"
    exit 1
fi

PARTITION_PREFIXES="0 1 2 3 4 5 6 7 8 9 a b c d e f"
S3_SYNC_CMD = "/usr/local/bin/aws s3 sync --copy-props metadata-directive --no-progress"
# TODO: Remove echo statements once we're sure this is working
echo parallel --line-buffer --jobs 16 $S3_SYNC_CMD s3://corpora-data-prod/ s3://corpora-data-${DEPLOYMENT_STAGE}/ --exclude "'*'" --include "'{}*'" --exclude "'*loom'" ::: $PARTITION_PREFIXES
echo parallel --line-buffer --jobs 16 $S3_SYNC_CMD s3://hosted-cellxgene-prod/ s3://hosted-cellxgene-${DEPLOYMENT_STAGE}/ --exclude "'*'" --include "'{}*'" ::: $PARTITION_PREFIXES
