#!/usr/bin/env bash

# Loads a pg_dump-produced data file into a PostgreSQL
# database. Updates table columns that contain deployment
# environment-specific URLs to point at the current environment.
#
# THIS IS DESTRUCTIVE! (It will not run in prod)
#
# Pre-requisities:
# - ssh tunnel to bastion is up
# - DEPLOYMENT_STAGE env var defined (must be 'dev' or 'staging')
#
# Usage:
# db_load.sh <pg_dump_file>
#
# where <pg_dump_file> is in "custom" format (pg_dump -Fc option); see backend/Makefile db/dump target


if [[ ! "${DEPLOYMENT_STAGE}" =~ ^(dev|staging) ]]; then
		echo ERROR: invalid value for destination deployment env: \"${DEPLOYMENT_STAGE}\";
		exit 1;
fi

DB_PW=`aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'`

INFILE=${1:?}
# TODO: Remove the suffix when ready for the real thing!
DB_NAME="corpora_${DEPLOYMENT_STAGE}_sync_test"
DB_USER=corpora_${DEPLOYMENT_STAGE}

PGPASSWORD=${DB_PW} pg_restore --clean --if-exists --no-owner --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} ${INFILE}

PGPASSWORD=${DB_PW} psql --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} -c "UPDATE dataset SET explorer_url = regexp_replace(explorer_url, '(https:\\/\\/)(.+?)(\\/.+)', '\\1cellxgene.${DEPLOYMENT_STAGE}.single-cell.czi.technology\\3') WHERE explorer_url IS NOT NULL"

PGPASSWORD=${DB_PW} psql --dbname=${DB_NAME} --host 0.0.0.0 --username ${DB_USER} -c "UPDATE dataset_artifact SET s3_uri = regexp_replace(s3_uri, '(s3:\\/\\/)([[:alpha:]]+-[[:alpha:]]+-)([[:alpha:]]+)(\\/.+)', '\\1\\2${DEPLOYMENT_STAGE}\\4') WHERE s3_uri IS NOT NULL"
