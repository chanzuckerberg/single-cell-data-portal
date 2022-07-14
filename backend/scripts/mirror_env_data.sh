#!/usr/bin/env bash

# Mirror S3 and RDS (postgres db) data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh
# Note: we run RDS mirror first to avoid possibility of referencing S3 data artifacts that are created after RDS dump
$SCRIPTS_DIR/mirror_rds_data.sh $@
$SCRIPTS_DIR/mirror_s3_data.sh $@


