#!/bin/bash
echo
echo " ====="
echo "| starting backend container"
echo " ====="
echo

# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  # Use locally-generated cert for HTTPS in containerized local environment; deployed envs use ELB
  HTTPS_CERT_AND_KEY="--certfile /tmp/pkcs12/server.crt --keyfile /tmp/pkcs12/server.key"
fi
# Note: Using just 1 worker for dev/test env. Multiple workers are used in deployment envs, as defined in Terraform code.
# Note: keep-alive timeout should always be greater than the idle timeout of the load balancer (60 seconds)
exec gunicorn ${HTTPS_CERT_AND_KEY} --worker-class gevent --workers 1 --bind 0.0.0.0:5000 backend.api_server.app:app \
  --max-requests 10000 --timeout 540 --keep-alive 61 --log-level info
