#!/bin/bash
echo
echo " ====="
echo "| starting backend container"
echo " ====="
echo

# Download WMG data snapshot to the local disk of the compute node on AWS
# This is done as optimization because retrieving data from local disk is
# significantly faster than retrieving data from S3
if [[ "${DEPLOYMENT_STAGE}" == "rdev" && -n "${REMOTE_DEV_PREFIX}" ]]; then
  echo "| Downloading WMG data snapshot for RDEV stack: ${REMOTE_DEV_PREFIX} from S3 to local disk..."
  strip_slash_remote_dev_prefix="${REMOTE_DEV_PREFIX//\//}" # strips ALL "/"
  echo "WMG_DISK_CACHE_PATH env var value is: ${WMG_DISK_CACHE_PATH}"
  echo aws s3 sync "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots" "/tmp/wmg_disk_cache/${strip_slash_remote_dev_prefix}/snapshots"
  aws s3 sync "s3://env-rdev-wmg/${strip_slash_remote_dev_prefix}/snapshots" "/tmp/wmg_disk_cache/${strip_slash_remote_dev_prefix}/snapshots"
elif [[ "${DEPLOYMENT_STAGE}" == "dev" || "${DEPLOYMENT_STAGE}" == "staging" || "${DEPLOYMENT_STAGE}" == "prod" ]]; then
  echo "| Downloading WMG data snapshot for deployment env: ${DEPLOYMENT_STAGE} from S3 to local disk..."
  echo "WMG_DISK_CACHE_PATH env var value is: ${WMG_DISK_CACHE_PATH}"
  echo aws s3 sync "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots" /tmp/wmg_disk_cache/snapshots
  aws s3 sync "s3://cellxgene-wmg-${DEPLOYMENT_STAGE}/snapshots" /tmp/wmg_disk_cache/snapshots
else
  echo "| Skipping downloading WMG data snapshot for deployment env: ${DEPLOYMENT_STAGE}..."
fi

echo "| Finished downloading WMG data snapshot from S3 for valid deployment environments"

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

