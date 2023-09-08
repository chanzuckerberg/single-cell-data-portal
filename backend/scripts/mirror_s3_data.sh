#!/usr/bin/env bash

# Mirror S3 data from a source deployment environment (usually production)
# to a specified destination deployment environment (dev or staging)
#
# This will *add* S3 objects to the dest env, but will not remove existing objects (see NOTE, below).
# The src S3 buckets will never be modified.

set -e

SCRIPTS_DIR=`dirname $0`
. $SCRIPTS_DIR/set_src_dest_envs.sh

# Set AWS_PROFILE according to SRC_ENV
if [[ $SRC_ENV == 'staging' ]]; then
  src_aws_profile=single-cell-dev
else
  src_aws_profile=single-cell-${SRC_ENV}
fi


# Copy WMB cube
echo "Copying WMG cube snapshot..."

# Note: "--copy-props metadata-directive" copies s3 object metadata, but not tags.
# We have not granted s3:GetObjectTagging perm to dev AWS account on prod buckets,
# so this avoids errors. As none of the s3 objects contain tags, this is acceptable.
AWS_OPTIONS="--copy-props metadata-directive --no-progress"
S3_COPY_CMD="aws s3 cp $AWS_OPTIONS"
S3_SYNC_CMD="aws s3 sync $AWS_OPTIONS"

latest_snapshot=$(aws s3 cp s3://cellxgene-wmg-${SRC_ENV}/latest_snapshot_run -)
if [[ $DEST_ENV == 'rdev' ]]; then
  s3_destination="env-rdev-wmg/${STACK}"
else
  s3_destination="cellxgene-wmg-${DEST_ENV}"
fi
echo -n $latest_snapshot | aws s3 cp - s3://${s3_destination}/latest_snapshot_run
$S3_SYNC_CMD s3://cellxgene-wmg-${SRC_ENV}/${latest_snapshot} s3://${s3_destination}/${latest_snapshot}


echo Mirroring S3 Dataset data from $SRC_ENV to $DEST_ENV...

if [[ $DEST_ENV == 'rdev' ]]; then

  # Clean up semaphores state directory ahead of using sem in case it was not cleared properly before
  rm -r ~/.parallel/semaphores/
  sample_data_name="schema_3_0_0"
  echo "Copying sample data for S3 Dataset assets for $(tr ',' '\n' <<< $COLLECTIONS | wc -l) Collections..."

  # Copy artifacts
  export AWS_PROFILE=$src_aws_profile  # For SRC_ENV
  DB_PW=$(aws secretsmanager get-secret-value --secret-id corpora/backend/${SRC_ENV}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string')
  # Select distinct local "directory" paths for assets, e.g. 's3://bucket/uuid/local.h5ad' becomes 's3://bucket/uuid/'
  uri_query="select distinct regexp_replace(uri, '(^.*\\/).*', '\\1') from \"DatasetArtifact\" where id in (select unnest(artifacts) from \"DatasetVersion\" where id in (select unnest(datasets) from \"CollectionVersion\" where collection_id in ('$(sed "s/,/','/g" <<< $COLLECTIONS)')))"
  DEPLOYMENT_STAGE=$SRC_ENV make db/tunnel/up
  uris=($(PGOPTIONS='-csearch_path=persistence_schema' PGPASSWORD=${DB_PW} psql --dbname corpora_${SRC_ENV} --username corpora_${SRC_ENV} --host 0.0.0.0 --csv --tuples-only -c "$uri_query"))
  export AWS_PROFILE=single-cell-dev  # For DEST_ENV

  for uri in "${uris[@]}"; do
    bucket=$(sed -E 's/s3:\/\/([^\/]+).*/\1/' <<< $uri)
    if [[ -n `grep 'hosted-cellxgene' <<< $bucket` ]]; then
      rdev_bucket_suffix="cellxgene"
    elif [[ -n `grep 'corpora-data' <<< $bucket` ]]; then
      rdev_bucket_suffix="artifacts"
    fi
    key_dir=$(sed -E 's/s3:\/\/([^\/]+)\/(.*)/\2/' <<< $uri)
    sem -j+0 $S3_SYNC_CMD s3://env-rdev-dataportal/s3/${rdev_bucket_suffix}/${sample_data_name}/ s3://env-rdev-${rdev_bucket_suffix}/${STACK}/${key_dir}
  done


  # Copy public assets
  export AWS_PROFILE=$src_aws_profile  # For SRC_ENV
  dv_id_query="select distinct unnest(datasets) from \"CollectionVersion\" where collection_id in ('$(sed "s/,/','/g" <<< $COLLECTIONS)')"
  dataset_version_ids=($(PGOPTIONS='-csearch_path=persistence_schema' PGPASSWORD=${DB_PW} psql --dbname corpora_${SRC_ENV} --username corpora_${SRC_ENV} --host 0.0.0.0 --csv --tuples-only -c "$dv_id_query"))
  DEPLOYMENT_STAGE=$SRC_ENV make db/tunnel/down  # db access no longer needed
  export AWS_PROFILE=single-cell-dev  # For DEST_ENV

  exts=("rds" "h5ad")
  for dv_id in "${dataset_version_ids[@]}"; do
    for ext in "${exts[@]}"; do
      uri="s3://dataset-assets-public-${SRC_ENV}/${dv_id}.${ext}"
      # Only copy if file exists in SRC_ENV public assets bucket (necessary to replicate tombstoned Dataset behavior)
      sem -j+0 "aws s3 ls $uri && $S3_COPY_CMD s3://env-rdev-dataportal/s3/datasets/${sample_data_name}.${ext} s3://env-rdev-datasets/${STACK}/${dv_id}.${ext}"
    done
  done
  sem --wait
else  # i.e., if DEST_ENV != 'rdev'
  # TODO: Add --delete once we confirm that is no data in the folders
  # that needs to be kept around (buckets are versioned, so we're not in
  # jeopardy of losing anything permanently). This would make the
  # operation more destructive, of course!

  PARALLEL_CMD="parallel --ungroup --jobs 16"  # '--ungroup' permits stdout to stream continuously from all jobs

  hex_chars="0 1 2 3 4 5 6 7 8 9 a b c d e f"  # Permits running 16 concurrent aws processes -- each one syncs only uuids starting with its hex char
  set -x
  echo -e "\n\n*** The following parallel aws s3 sync commands may appear to hang, but they are working. Patience is a virtue. ***\n\n"
  $PARALLEL_CMD $S3_SYNC_CMD s3://corpora-data-${SRC_ENV}/ s3://corpora-data-${DEST_ENV}/ --exclude "'*'" --include "'{}*'" --exclude "'*loom'" ::: $hex_chars
  $PARALLEL_CMD $S3_SYNC_CMD s3://hosted-cellxgene-${SRC_ENV}/ s3://hosted-cellxgene-${DEST_ENV}/ --exclude "'*'" --include "'{}*'" ::: $hex_chars
  # We need to assume the sync-datasets-{dev|staging} role, which has the necessary permissions to perform a s3 sync from
  # {dev|staging|prod} to {dev|staging} **for the public-access datasets bucket**.
  ACCOUNT_ID=`aws sts get-caller-identity | jq -r '.Arn' | grep -o '[0-9]*'`  # First get AWS account id for dev account
  SYNC_ROLE_CREDENTIALS=`aws sts assume-role --role-arn arn:aws:iam::${ACCOUNT_ID}:role/sync-datasets-dev \
                                             --role-session-name sync-datasets \
                         | jq -r '.Credentials'`
  # Execute sync for datasets bucket with sync role credentials
  AWS_SESSION_TOKEN=`jq -r '.SessionToken' <<< $SYNC_ROLE_CREDENTIALS` \
  AWS_SECRET_ACCESS_KEY=`jq -r '.SecretAccessKey' <<< $SYNC_ROLE_CREDENTIALS` \
  AWS_ACCESS_KEY_ID=`jq -r '.AccessKeyId' <<< $SYNC_ROLE_CREDENTIALS` \
  $PARALLEL_CMD $S3_SYNC_CMD s3://dataset-assets-public-${SRC_ENV}/ s3://dataset-assets-public-${DEST_ENV}/ --exclude "'*'" --include "'{}*'" ::: $hex_chars
  set +x
fi
