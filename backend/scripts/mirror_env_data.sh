#!/usr/bin/env bash

# Mirror S3 and RDS (postgres db) data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)
#
# THIS IS DESTRUCTIVE for the destination env! The source env will
# never be modified, but the dest env's data will be replaced.

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh
export NO_PROMPT=1
# Note: Run RDS mirroring before S3 mirroring. This avoids edge case
# of having dangling references to S3 data artifacts in the dest env
# db; this can occur if new artifacts are created and recorded in the
# src env after the S3 mirroring object set has been computed.
#$SCRIPTS_DIR/mirror_rds_data.sh $@
if [[ ! $DEST_ENV == 'rdev' ]]; then
#  $SCRIPTS_DIR/mirror_s3_data.sh $@
  echo "ok"
elif [[ ! -z "$COLLECTIONS" ]]; then
  echo -e "\nCopying assets for the following Collections:\n$(tr ',' '\n' <<< $COLLECTIONS)"
  $SCRIPTS_DIR/mirror_s3_data.sh $@
else
  echo -e "\nDEST_ENV is set to rdev -- will NOT copy s3 assets"
fi

echo "Done!"
