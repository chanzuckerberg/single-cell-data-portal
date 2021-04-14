#!/bin/bash
export AWS_REGION=us-west-2
export AWS_DEFAULT_REGION=us-west-2
export AWS_ACCESS_KEY_ID=nonce
export AWS_SECRET_ACCESS_KEY=nonce

export FRONTEND_URL=http://frontend.corporanet.local:3000
export BACKEND_URL=http://backend.corporanet.local:5000

# NOTE: This script is intended to run INSIDE the dockerized dev environment!
# If you need to run it directly on your laptop for some reason, change
# localstack below to localhost
export LOCALSTACK_URL=http://localstack.corporanet.local:4566
# How the backend can reach the OIDC idp
export OIDC_INTERNAL_URL=http://oidc.corporanet.local
# How a web browser can reach the OIDC idp
export OIDC_BROWSER_URL=https://oidc.corporanet.local:8443

# Get test credentials from oauth/user.json
oauth_file="./users.json"
oauth_user=$(cat ${oauth_file} | jq '.[0]')

export TEST_ACCOUNT_USERNAME=$(jq '.Username' <<< "${oauth_user}")
export TEST_ACCOUNT_PASSWORD=$(jq '.Password' <<< "${oauth_user}")

echo -n "waiting for localstack to be ready: "
until $(curl --output /dev/null --silent --head ${LOCALSTACK_URL}); do
    echo -n '.'
    sleep 1
done
echo " done"

echo "Creating secretsmanager secrets"
local_aws="aws --endpoint-url=${LOCALSTACK_URL}"
${local_aws} s3api create-bucket --bucket corpora-data-dev &>/dev/null || true
${local_aws} s3api create-bucket --bucket artifact-bucket &>/dev/null || true
${local_aws} s3api create-bucket --bucket cellxgene-bucket &>/dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/test/auth0-secret &>/dev/null || true
${local_aws} secretsmanager create-secret --name corpora/cicd/test/auth0-secret &>/dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/test/database_local &>/dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/test/config &>/dev/null || true

echo "Creating default state machine"
${local_aws} iam create-role --role-name StepRole --assume-role-policy-document '{"Version": "2012-10-17", "Statement": [ { "Effect": "Allow", "Principal": { "Service": "states.amazonaws.com" }, "Action": "sts:AssumeRole" }' || true
${local_aws} stepfunctions create-state-machine --name uploader-dev-sfn --definition '{"StartAt": "noop", "States": {"noop": {"Type": "Pass", "Result": {}, "ResultPath": "$.coords", "End": true}}}' --role-arn arn:aws:iam::000000000000:role/StepRole || true

echo "Updating secrets"
${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/auth0-secret --secret-string '{
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
    "redirect_to_frontend": "'"${FRONTEND_URL}"'",
    "test_account_username": '"${TEST_ACCOUNT_USERNAME}"',
    "test_account_password": '"${TEST_ACCOUNT_PASSWORD}"'
}' || true


# TODO: python3 -m unittest tests.unit.backend.corpora.common.test_authorizer.TestAuthorizer.test_invalid_token
${local_aws} secretsmanager update-secret --secret-id corpora/cicd/test/auth0-secret --secret-string '{
    "client_id": "",
    "client_secret": "",
    "audience": "",
    "grant_type": ""
}' || true

${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/database_local --secret-string '{"database_uri": "postgresql://corpora:test_pw@database.corporanet.local:5432"}' || true
${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/config --secret-string '{"upload_sfn_arn": "arn:aws:states:us-west-2:000000000000:stateMachine:uploader-dev-sfn"}' || true

# Make a 1mb data file
echo "Writing test file to s3"
dd if=/dev/zero of=fake-h5ad-file.h5ad bs=1024 count=1024 &>/dev/null
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
