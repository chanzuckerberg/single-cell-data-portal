#!/bin/bash
echo
echo " ====="
echo "| starting backend container"
echo " ====="
echo

# Download WMG data snapshot to a mounted filesystem of the compute node on AWS
# This is done as optimization because retrieving data from local disk is
# significantly faster than retrieving data from S3
WMG_SNAPSHOT_FS_CACHE_ROOT_PATH="/single-cell-data-portal/wmg_snapshot_cache"
LATEST_READER_SNAPSHOT_SCHEMA_VERSION="v3"

echo "| ENV VAR DOWNLOAD_WMG_DATA_TO_DISK: ${DOWNLOAD_WMG_DATA_TO_DISK}"

if [[ "${DOWNLOAD_WMG_DATA_TO_DISK}" == "false" ]]; then
  echo "| Skipping downloading WMG data snapshot because DOWNLOAD_WMG_DATA_TO_DISK is set to false"
elif [[ "${DEPLOYMENT_STAGE}" == "rdev" && -n "${REMOTE_DEV_PREFIX}" ]]; then
  echo "| Downloading WMG data snapshot for RDEV stack: ${REMOTE_DEV_PREFIX} from S3 to filesystem path: ${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}"

  strip_slash_remote_dev_prefix="${REMOTE_DEV_PREFIX//\//}" # strips ALL "/"

  latest_snapshot_identifier=$(aws s3 cp "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" -)

  echo aws s3 sync "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}"
  aws s3 sync "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}"

  echo aws s3 cp "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier"
  aws s3 cp "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/${strip_slash_remote_dev_prefix}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier"

elif [[ "${DEPLOYMENT_STAGE}" == "dev" || "${DEPLOYMENT_STAGE}" == "staging" || "${DEPLOYMENT_STAGE}" == "prod" ]]; then
  echo "| Downloading WMG data snapshot for deployment env: ${DEPLOYMENT_STAGE} from S3 to filesystem path: ${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}"

  latest_snapshot_identifier=$(aws s3 cp "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" -)

  echo aws s3 sync "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}"
  aws s3 sync "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/${latest_snapshot_identifier}"

  echo aws s3 cp "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier"
  aws s3 cp "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier" "${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH}/snapshots/${LATEST_READER_SNAPSHOT_SCHEMA_VERSION}/latest_snapshot_identifier"

else
  echo "| Skipping downloading WMG data snapshot for deployment env: ${DEPLOYMENT_STAGE}..."
fi

echo "| Finished downloading WMG data snapshot from S3 to filesystem path: ${WMG_SNAPSHOT_FS_CACHE_ROOT_PATH} for valid deployment environments"

# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  # Use locally-generated cert for HTTPS in containerized local environment; deployed envs use ELB
  HTTPS_CERT_AND_KEY="--certfile /tmp/pkcs12/server.crt --keyfile /tmp/pkcs12/server.key"
fi

# TODO: Fix problem with invoking gunicorn under ddtrace-run. Also make sure when
# DEPLOYMENT_STAGE=test, gunicorn is not invoked under ddtrace-run.
# See ticket: https://github.com/chanzuckerberg/single-cell-data-portal/issues/5819

export DD_GEVENT_PATCH_ALL=true
# Note: Using just 1 worker for dev/test env. Multiple workers are used in deployment envs, as defined in Terraform code.
# Note: keep-alive timeout should always be greater than the idle timeout of the load balancer (60 seconds)

echo "starting gunicorn server"
exec gunicorn ${HTTPS_CERT_AND_KEY} --worker-class gevent --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app \
  --max-requests 10000 --timeout 180 --keep-alive 61 --log-level info
