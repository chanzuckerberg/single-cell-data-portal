#!/bin/bash
export AWS_REGION=us-west-2
export AWS_DEFAULT_REGION=us-west-2
export AWS_ACCESS_KEY_ID=nonce
export AWS_SECRET_ACCESS_KEY=nonce

export FRONTEND_URL=https://frontend.corporanet.local:3000
export BACKEND_URL=https://backend.corporanet.local:5000

# NOTE: This script is intended to run INSIDE the dockerized dev environment!
# If you need to run it directly on your laptop for some reason, change
# localstack below to localhost
export LOCALSTACK_URL=http://localstack.corporanet.local:4566
# How the backend can reach the OIDC idp
export OIDC_INTERNAL_URL=http://auth.corporanet.local
# How a web browser can reach the OIDC idp
export OIDC_BROWSER_URL=https://auth.corporanet.local

# Get test credentials from oauth/users.json
oauth_file="oauth/users.json"
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
${local_aws} secretsmanager create-secret --name corpora/backend/test/database &>/dev/null || true
${local_aws} secretsmanager create-secret --name corpora/backend/test/config &>/dev/null || true

echo "Creating default state machine"
${local_aws} iam create-role --role-name StepRole --assume-role-policy-document '{"Version": "2012-10-17", "Statement": [ { "Effect": "Allow", "Principal": { "Service": "states.amazonaws.com" }, "Action": "sts:AssumeRole" }' || true
${local_aws} stepfunctions create-state-machine --name uploader-dev-sfn --definition '{"StartAt": "noop", "States": {"noop": {"Type": "Pass", "Result": {}, "ResultPath": "$.coords", "End": true}}}' --role-arn arn:aws:iam::000000000000:role/StepRole || true

echo "Updating secrets"
${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/auth0-secret --secret-string '{
    "client_id": "local-client-id",
    "client_secret": "local-client-secret",
    "audience": "local-client-id",
    "api_audience": "local-client-id",
    "code_challenge_method": "S256",
    "api_authorize_url": "'"${OIDC_BROWSER_URL}"'/connect/authorize",
    "api_base_url": "'"${OIDC_BROWSER_URL}"'",
    "api_token_url": "'"${OIDC_INTERNAL_URL}"'/connect/token",
    "api_userinfo_url": "'"${OIDC_INTERNAL_URL}"'/connect/userinfo",
    "internal_url": "'"${OIDC_INTERNAL_URL}"'",
    "cookie_name": "cxguser",
    "callback_base_url": "'"${BACKEND_URL}"'",
    "redirect_to_frontend": "'"${FRONTEND_URL}"'",
    "test_account_username": '"${TEST_ACCOUNT_USERNAME}"',
    "test_account_password": '"${TEST_ACCOUNT_PASSWORD}"',
    "domain":"https://localhost:1234",
    "api_key_secret": "a random secret",
    "days_to_live": 60,
    "mgmt_client_id": "test_mgmt_client_id",
    "mgmt_client_secret": "test_mgmt_client_secret",
    "api_key_connection_name": "api-key-database",
    "auth0_domain" :"localhost",
    "curation_audience": "localhost/curation"
}' || true


# TODO: python3 -m unittest tests.unit.backend.common.test_authorizer.TestAuthorizer.test_invalid_token
${local_aws} secretsmanager update-secret --secret-id corpora/cicd/test/auth0-secret --secret-string '{
    "client_id": "",
    "client_secret": "",
    "audience": "",
    "grant_type": ""
}' || true

${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/database --secret-string '{"database_uri":
 "postgresql://corpora:test_pw@database.corporanet.local:5432"}' || true
${local_aws} secretsmanager update-secret --secret-id corpora/backend/test/config --secret-string '{"upload_sfn_arn": "arn:aws:states:us-west-2:000000000000:stateMachine:uploader-dev-sfn", "curator_role_arn":"test_curation_role"}' || true

# Make a 1mb data file
echo "Writing test file to s3"
dd if=/dev/zero of=fake-h5ad-file.h5ad bs=1024 count=1024 &>/dev/null
${local_aws} s3 cp fake-h5ad-file.h5ad s3://corpora-data-dev/
rm fake-h5ad-file.h5ad

echo "Populating test db"
export BOTO_ENDPOINT_URL=${LOCALSTACK_URL}
cd $(dirname ${BASH_SOURCE[0]})/..
python3 -m scripts.populate_db

# Make a WMG snapshot
echo "Setting up WMG snapshot data"
tmp_snapshot_dir=`mktemp -d`
snapshot_identifier='dummy-snapshot'
wmg_bucket="wmg-test"
wmg_config_secret_name="corpora/backend/test/wmg_config"
python3 -m tests.unit.backend.wmg.fixtures.test_snapshot ${tmp_snapshot_dir}
${local_aws} s3api create-bucket --bucket ${wmg_bucket} &>/dev/null || true
${local_aws} s3 sync --delete --quiet ${tmp_snapshot_dir} s3://${wmg_bucket}/$snapshot_identifier/
echo $snapshot_identifier | ${local_aws} s3 cp --quiet - s3://${wmg_bucket}/latest_snapshot_identifier
${local_aws} secretsmanager create-secret --name ${wmg_config_secret_name} &>/dev/null || true
${local_aws} secretsmanager update-secret --secret-id ${wmg_config_secret_name} --secret-string "{\"bucket\": \"${wmg_bucket}\"}" || true

# Access NCBI API key
echo "Creating NCBI API key secret"
ncbi_api_key="test-key"
gene_info_config_secret_name="corpora/backend/test/gene_info_config"
${local_aws} secretsmanager create-secret --name ${gene_info_config_secret_name} &>/dev/null || true
${local_aws} secretsmanager update-secret --secret-id ${gene_info_config_secret_name} --secret-string "{\"ncbi_api_key\": \"${ncbi_api_key}\"}" || true

echo
echo "Dev env is up and running!"
echo "  Frontend: ${FRONTEND_URL}"
echo "  Backend: ${BACKEND_URL}"
