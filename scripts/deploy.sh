#!/usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage:"
  echo "  scripts/deploy.sh rdev --pr <number> [--image-tag <tag>]"
  echo "  scripts/deploy.sh staging [--image-tag <tag>]"
  echo "  scripts/deploy.sh prod [--image-tag <tag>]"
}

fail() {
  echo "$1" >&2
  exit 1
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || fail "$1 is required"
}

clear_aws_credentials() {
  unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN
}

assume_role() {
  local account_id="$1"
  clear_aws_credentials
  credentials=$(AWS_PROFILE=czi-id aws sts assume-role \
    --role-arn "arn:aws:iam::${account_id}:role/tfe-si" \
    --role-session-name single-cell-data-portal-deploy)
  export AWS_ACCESS_KEY_ID
  AWS_ACCESS_KEY_ID=$(jq -r '.Credentials.AccessKeyId' <<<"$credentials")
  export AWS_SECRET_ACCESS_KEY
  AWS_SECRET_ACCESS_KEY=$(jq -r '.Credentials.SecretAccessKey' <<<"$credentials")
  export AWS_SESSION_TOKEN
  AWS_SESSION_TOKEN=$(jq -r '.Credentials.SessionToken' <<<"$credentials")
}

assume_target_role() {
  assume_role "$target_account"
}

images_exist() {
  local image_tag="$1"
  local missing=0
  local repositories=(
    corpora-frontend
    corpora-backend
    corpora-backend-de
    corpora-backend-wmg
    corpora-upload-failures
    corpora-upload-success
    corpora-upload
    wmg-processing
    cellguide-pipeline
  )

  assume_role 699936264352
  for repository in "${repositories[@]}"; do
    if ! aws ecr describe-images \
      --repository-name "$repository" \
      --image-ids "imageTag=$image_tag" >/dev/null 2>&1; then
      missing=1
    fi
  done
  clear_aws_credentials
  export AWS_PROFILE=czi-id

  return "$missing"
}

wait_for_image_build() {
  local deployment_id="$1"
  local run_id=""

  for _ in {1..60}; do
    runs=$(gh run list \
      --workflow "$workflow" \
      --branch "$source_ref" \
      --event workflow_dispatch \
      --limit 20 \
      --json databaseId,displayTitle)
    run_id=$(jq -r \
      --arg deployment_id "$deployment_id" \
      'map(select(.displayTitle == $deployment_id))[0].databaseId // empty' \
      <<<"$runs")
    if [[ -n "$run_id" ]]; then
      echo "$run_id"
      return
    fi
    sleep 2
  done

  fail "The image build workflow run did not appear"
}

run_database_migration() {
  task_definition=$(terraform output -raw migrate_db_task_definition_arn)
  assume_target_role

  config=$(aws secretsmanager get-secret-value \
    --secret-id "$config_secret" \
    --query SecretString \
    --output text)
  cluster=$(jq -r '.cluster_arn' <<<"$config")
  subnets=$(jq -c '.private_subnets' <<<"$config")
  security_groups=$(jq -c '.security_groups' <<<"$config")

  run_result=$(aws ecs run-task \
    --cluster "$cluster" \
    --task-definition "$task_definition" \
    --launch-type FARGATE \
    --network-configuration \
      "awsvpcConfiguration={subnets=${subnets},securityGroups=${security_groups},assignPublicIp=DISABLED}")
  failures=$(jq '.failures | length' <<<"$run_result")
  [[ "$failures" == "0" ]] || fail "The database migration task did not start"
  task_arn=$(jq -r '.tasks[0].taskArn' <<<"$run_result")
  [[ "$task_arn" != "null" ]] || fail "The database migration task did not return an ARN"

  while true; do
    task_status=$(aws ecs describe-tasks \
      --cluster "$cluster" \
      --tasks "$task_arn" \
      --query 'tasks[0].lastStatus' \
      --output text)
    [[ "$task_status" != "None" ]] || fail "The database migration task could not be found"
    [[ "$task_status" != "STOPPED" ]] || break
    sleep 15
  done

  exit_codes=$(aws ecs describe-tasks \
    --cluster "$cluster" \
    --tasks "$task_arn" \
    --query 'tasks[0].containers[].exitCode' \
    --output text)
  [[ -n "$exit_codes" && "$exit_codes" != "None" ]] || fail "The database migration task did not return an exit code"
  for exit_code in $exit_codes; do
    [[ "$exit_code" == "0" ]] || fail "The database migration task failed with exit code $exit_code"
  done

  clear_aws_credentials
  export AWS_PROFILE=czi-id
}

invalidate_cloudfront() {
  [[ "$target" != "rdev" ]] || return

  assume_target_role
  distributions=$(aws cloudfront list-distributions)
  distribution_ids=$(jq -c \
    --arg domain_name "$domain_name" \
    --arg domain_alias "$domain_alias" \
    '[.DistributionList.Items[]? | select(any(.Origins.Items[]?; .DomainName | contains($domain_name))) | select((.Aliases.Items // []) | index($domain_alias)) | .Id]' \
    <<<"$distributions")
  distribution_count=$(jq 'length' <<<"$distribution_ids")
  [[ "$distribution_count" == "1" ]] || fail "Expected one matching CloudFront distribution, found $distribution_count"
  distribution_id=$(jq -r '.[0]' <<<"$distribution_ids")
  aws cloudfront create-invalidation \
    --distribution-id "$distribution_id" \
    --paths /index.html

  clear_aws_credentials
  export AWS_PROFILE=czi-id
}

target="${1:-}"
[[ -n "$target" ]] || {
  usage
  exit 1
}
shift

pr_number=""
image_tag_override=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pr)
      [[ -n "${2:-}" ]] || {
        usage
        exit 1
      }
      pr_number="$2"
      shift 2
      ;;
    --image-tag)
      [[ -n "${2:-}" ]] || {
        usage
        exit 1
      }
      image_tag_override="$2"
      shift 2
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

if [[ "$target" == "rdev" ]]; then
  [[ -n "$pr_number" ]] || {
    usage
    exit 1
  }
elif [[ "$target" != "staging" && "$target" != "prod" ]]; then
  usage
  exit 1
elif [[ -n "$pr_number" ]]; then
  usage
  exit 1
fi
[[ -z "$image_tag_override" || "$image_tag_override" =~ ^sha-[0-9a-f]{8}$ ]] || fail "Image tags must use sha- followed by eight lowercase hexadecimal characters"

for command in aws gh git jq terraform; do
  require_command "$command"
done

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
[[ -z "$(git status --porcelain)" ]] || fail "The working tree must be clean"

repo=$(gh repo view --json nameWithOwner --jq '.nameWithOwner')
[[ "$repo" == "chanzuckerberg/single-cell-data-portal" ]] || fail "Run this script from chanzuckerberg/single-cell-data-portal"

git fetch origin

case "$target" in
  rdev)
    pr=$(gh pr view "$pr_number" \
      --json headRefName,headRefOid,isCrossRepository,state)
    [[ "$(jq -r '.state' <<<"$pr")" == "OPEN" ]] || fail "Pull request $pr_number is not open"
    [[ "$(jq -r '.isCrossRepository' <<<"$pr")" == "false" ]] || fail "Cross-repository pull requests are not supported"
    source_ref=$(jq -r '.headRefName' <<<"$pr")
    source_sha=$(jq -r '.headRefOid' <<<"$pr")
    terraform_root="$repo_root/.happy/terraform/envs/rdev"
    target_account=699936264352
    config_secret=happy/env-rdev-config
    ;;
  staging)
    source_ref=main
    source_sha=$(git rev-parse origin/main)
    terraform_root="$repo_root/.happy/terraform/envs/stage"
    target_account=699936264352
    config_secret=happy/env-stage-config
    domain_name=frontend.stage.single-cell.czi.technology
    domain_alias=cellxgene.staging.single-cell.czi.technology
    ;;
  prod)
    source_ref=main
    source_sha=$(git rev-parse origin/main)
    terraform_root="$repo_root/.happy/terraform/envs/prod"
    target_account=231426846575
    config_secret=happy/env-prod-config
    domain_name=frontend.production.single-cell.czi.technology
    domain_alias=cellxgene.cziscience.com
    ;;
esac

[[ "$(git rev-parse HEAD)" == "$source_sha" ]] || fail "Check out commit $source_sha before deploying $target"

terraform_version=$(terraform version -json | jq -r '.terraform_version')
[[ "$terraform_version" == 1.3.* ]] || fail "Terraform 1.3 is required"

export AWS_PROFILE=czi-id
export AWS_REGION=us-west-2
clear_aws_credentials
aws sts get-caller-identity >/dev/null

workflow=build-images-and-create-deployment.yml
image_tag="${image_tag_override:-sha-${source_sha:0:8}}"
if images_exist "$image_tag"; then
  echo "Reusing existing images tagged $image_tag"
elif [[ -n "$image_tag_override" ]]; then
  fail "Not all nine images exist with tag $image_tag"
else
  deployment_id="deploy-${target}-${source_sha:0:8}-$(date +%s)-$$"
  gh workflow run "$workflow" \
    --ref "$source_ref" \
    --field deployment_id="$deployment_id"
  run_id=$(wait_for_image_build "$deployment_id")
  gh run watch "$run_id" --exit-status
  image_sha=$(gh run view "$run_id" --json headSha --jq '.headSha')
  [[ "$image_sha" == "$source_sha" ]] || fail "The image build used $image_sha instead of $source_sha"
  image_tag="sha-${image_sha:0:8}"
fi

cd "$terraform_root"
terraform init -reconfigure
if [[ "$target" != "rdev" ]]; then
  state_count=$(terraform state list | wc -l | tr -d ' ')
  [[ "$state_count" != "0" ]] || fail "The $target state is empty"
fi

plan_file=$(mktemp "${TMPDIR:-/tmp}/single-cell-data-portal-${target}.XXXXXX")
trap 'rm -f "$plan_file"; clear_aws_credentials' EXIT
terraform plan -var "image_tag=$image_tag" -out "$plan_file"
terraform show "$plan_file"

read -r -p "Apply $image_tag to $target? [y/N] " answer
[[ "$answer" == "y" || "$answer" == "Y" ]] || {
  echo "Deployment cancelled"
  exit 0
}

terraform apply "$plan_file"
run_database_migration
invalidate_cloudfront
terraform output

echo "Deployed $image_tag to $target"
