#!/bin/bash
docker restart $(docker ps -f name=backend -q)

export LOCALSTACK_URL=http://localstack.corporanet.local:4566
local_aws="aws --endpoint-url=${LOCALSTACK_URL}"

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
${local_aws} secretsmanager update-secret --secret-id ${wmg_config_secret_name} --secret-string '{"bucket": "wmg-test", "tiledb_config_overrides": "{\"py.init_buffer_bytes\": 536870912, 
\"sm.tile_cache_size\": 134217728, \"sm.mem.total_budget\": 1073741824, \"sm.memory_budget\": 536870912, \"sm.memory_budget_var\": 1073741824}"}'

echo
echo "Backend restarted!"
