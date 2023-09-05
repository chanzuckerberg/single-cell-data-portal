#!/usr/bin/env bash

# Mirror S3 data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)
#
# This will *add* S3 objects to the dest env, but will not remove existing objects (see NOTE, below).
# The src S3 buckets will never be modified.

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh

DEPLOYMENT_STAGE=$SRC_ENV

echo Copying WMG cube snapshot...

latest_snapshot=$(aws s3 cp s3://cellxgene-wmg-${SRC_ENV}/latest_snapshot_run -)
if [[ $DEST_ENV == 'rdev' ]]; then
  s3_destination="env-rdev-wmg/${STACK}"
else
  s3_destination="cellxgene-wmg-${DEST_ENV}"
fi
echo -n $latest_snapshot | aws s3 cp - s3://${s3_destination}/latest_snapshot_run
aws s3 cp --recursive s3://cellxgene-wmg-${SRC_ENV}/${latest_snapshot} s3://${s3_destination}/${latest_snapshot}


echo Mirroring S3 Dataset data from $SRC_ENV to $DEST_ENV...

if [[ $DEST_ENV == 'rdev' ]]; then
  DEPLOYMENT_STAGE=$DEPLOYMENT_STAGE make db/tunnel/up
  echo "Copying S3 Dataset assets for ${#COLLECTIONS[@]} Collections..."
  DB_PW=$(aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string')

  # Copy artifacts
  uri_query="select uri from \"DatasetArtifact\" where id in (select unnest(artifacts) from \"DatasetVersion\" where id in (select unnest(datasets) from \"CollectionVersion\" where collection_id in ('$(sed "s/,/','/g" <<< $COLLECTIONS)')))"
  uris=($(PGOPTIONS='-csearch_path=persistence_schema' PGPASSWORD=${DB_PW} psql --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0 --csv --tuples-only -c "$uri_query"))
  for uri in "${uris[@]}"; do
    bucket=$(sed -E 's/s3:\/\/([^\/]+).*/\1/' <<< $uri)
    if [[ ! -z `grep 'hosted-cellxgene' <<< $bucket` ]]; then
      rdev_bucket_suffix="cellxgene"
    elif [[ ! -z `grep 'corpora-data' <<< $bucket` ]]; then
      rdev_bucket_suffix="artifacts"
    fi
    key=$(sed -E 's/s3:\/\/([^\/]+)\/(.*)/\2/' <<< $uri)
    recursive_flag=$([[ ! "${key: -1}" == "/" ]] || echo "--recursive ")
    aws s3 cp ${recursive_flag}${uri} s3://env-rdev-${rdev_bucket_suffix}/${STACK}/${key}
  done

  # Copy public assets
  dv_id_query="select unnest(datasets) from \"CollectionVersion\" where collection_id in ('$(sed "s/,/','/g" <<< $COLLECTIONS)')"
  dataset_version_ids=($(PGOPTIONS='-csearch_path=persistence_schema' PGPASSWORD=${DB_PW} psql --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0 --csv --tuples-only -c "$dv_id_query"))
  DEPLOYMENT_STAGE=$DEPLOYMENT_STAGE make db/tunnel/down
  exts=("rds" "h5ad")
  for dv_id in "${dataset_version_ids[@]}"; do
    for ext in "${exts[@]}"; do
      aws s3 cp s3://dataset-assets-public-${DEPLOYMENT_STAGE}/${dv_id}.${ext} s3://env-rdev-datasets/${STACK}/${dv_id}.${ext}
    done
  done

  exit 0
fi

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
