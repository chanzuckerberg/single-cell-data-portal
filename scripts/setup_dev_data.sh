#!/bin/bash
export AWS_REGION=us-west-2
export AWS_DEFAULT_REGION=us-west-2
export AWS_ACCESS_KEY_ID=nonce
export AWS_SECRET_ACCESS_KEY=nonce

export FRONTEND_URL=http://localhost:8000
export BACKEND_URL=http://localhost:5000

# NOTE: This script is intended to run INSIDE the dockerized dev environment!
# If you need to run it directly on your laptop for some reason, change
# localstack below to localhost
export LOCALSTACK_URL=http://localstack:4566
# How the backend can reach the OIDC idp
export OIDC_INTERNAL_URL=http://oidc
# How a web browser can reach the OIDC idp
export OIDC_BROWSER_URL=https://localhost:8443

echo -n "waiting for localstack to be ready: "
until $(curl --output /dev/null --silent --head ${LOCALSTACK_URL}); do
    echo -n '.'
    sleep 1
done
echo " done"

echo "Creating secretsmanager secrets"
local_aws="aws --endpoint-url=${LOCALSTACK_URL}"
${local_aws} s3api create-bucket --bucket corpora-data-dev &> /dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/dev/auth0-secret &> /dev/null || true
${local_aws} secretsmanager create-secret --name corpora/cicd/test/auth0-secret &> /dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/dev/database_local &> /dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/test/database_local &> /dev/null || true


echo "Updating secrets"
${local_aws} secretsmanager update-secret --secret-id corpora/backend/dev/auth0-secret --secret-string '{
    "client_id": "local-client-id",
    "client_secret": "local-client-secret",
    "audience": "local-client-id",
    "code_challenge_method": "S256",
    "api_authorize_url": "'"${OIDC_BROWSER_URL}"'/connect/authorize",
    "api_base_url": "'"${OIDC_BROWSER_URL}"'",
    "api_token_url": "'"${OIDC_INTERNAL_URL}"'/connect/token",
    "internal_url": "'"${OIDC_INTERNAL_URL}"'",
    "cookie_name": "cxguser",
    "callback_base_url": "'"${BACKEND_URL}"'",
    "redirect_to_frontend": "'"${FRONTEND_URL}"'"
}' || true

# TODO: python3 -m unittest tests.unit.backend.corpora.common.test_authorizer.TestAuthorizer.test_invalid_token
${local_aws} secretsmanager update-secret --secret-id corpora/cicd/test/auth0-secret --secret-string '{
    "client_id": "",
    "client_secret": "",
    "audience": "",
    "grant_type": ""
}' || true

${local_aws} secretsmanager update-secret --secret-id corpora/backend/dev/database_local --secret-string '{"database_uri": "postgresql://corpora:test_pw@database:5432"}' || true
${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/database_local --secret-string '{"database_uri": "postgresql://corpora:test_pw@database:5432"}' || true

# Make a 1mb data file
echo "Writing test file to s3"
dd if=/dev/zero of=fake-h5ad-file.h5ad bs=1024 count=1024 &> /dev/null
${local_aws} s3 cp fake-h5ad-file.h5ad s3://corpora-data-dev/
rm fake-h5ad-file.h5ad

echo "Populating test db"
export CORPORA_LOCAL_DEV=true
export BOTO_ENDPOINT_URL=${LOCALSTACK_URL}
python3 $(dirname ${BASH_SOURCE[0]})/populate_db.py
echo
echo "Dev env is up and running!"
echo "  Frontend: ${FRONTEND_URL}"
echo "  Backend: ${BACKEND_URL}"
