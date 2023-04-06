#!/usr/bin/env bash

# Migrate S3 data from corpora-data-<env> bucket to dataset-assets-public-<env> bucket
set -e
export ENV=${1}

echo mirroring data in $ENV
if [[ "$ENV" == dev || "$ENV" == staging ]]; then
  export AWS_PROFILE=single-cell-dev
elif [[ "$ENV" == prod ]]; then
  export AWS_PROFILE=single-cell-prod
else
  echo "invalid env input"
  exit 1
fi

aws s3 ls s3://corpora-data-${ENV}/ | awk '{print $2}' | sed 's:/*$::' > filenames.txt
for file in `cat filenames.txt`; do
  echo "dry run copying s3://corpora-data-${ENV}/$file/local.h5ad TO s3://dataset-assets-public-${ENV}/$file.h5ad"
  # aws s3 cp s3://corpora-data-${ENV}/$file/local.h5ad s3://dataset-assets-public-${ENV}/$file.h5ad
  # aws s3 cp s3://corpora-data-${ENV}/$file/local.rds s3://dataset-assets-public-${ENV}/$file.rds
done


