# Manual deployment

The data portal no longer deploys through Happy or Terraform Enterprise (TFE). A core infrastructure engineer must run Terraform from a local checkout for a break-glass deployment.

## Prerequisites

- Terraform 1.3.0
- AWS access through the `czi-id` profile
- permission to assume `tfe-si` in the target account and in `core-platform-prod`
- an image tag that exists in every required Elastic Container Registry (ECR) repository
- a clean checkout of the commit to deploy

The deployment roots are:

| Environment | Directory | AWS account | Stack |
| --- | --- | --- | --- |
| dev | `.happy/terraform/envs/dev` | `699936264352` | `devstack` |
| rdev | `.happy/terraform/envs/rdev` | `699936264352` | `rdevstack` |
| staging | `.happy/terraform/envs/stage` | `699936264352` | `stagestack` |
| prod | `.happy/terraform/envs/prod` | `231426846575` | `prodstack` |

The rdev root manages one static break-glass stack. It does not create a stack for each pull request. Its first apply creates the stack because no rdev state existed when the TFE states moved to Amazon Simple Storage Service (S3).

## Plan

Set the environment directory and immutable image tag:

```bash
export AWS_PROFILE=czi-id
export AWS_REGION=us-west-2
export TF_ROOT=.happy/terraform/envs/dev
export IMAGE_TAG=sha-01234567

cd "$TF_ROOT"
terraform init -reconfigure
terraform state list
terraform plan -var "image_tag=$IMAGE_TAG" -out deploy.tfplan
terraform show deploy.tfplan
```

For staging, use `stage` as the directory name. Stop if `terraform state list` is empty for dev, staging or prod. An empty state means Terraform is using the wrong state key and the plan will try to recreate the environment.

Review replacements and deletions before continuing. The dev state predates the current configuration and can contain substantial drift. Do not include unrelated drift in a break-glass deployment.

## Apply

Apply only the saved plan that you reviewed:

```bash
terraform apply deploy.tfplan
terraform output
rm deploy.tfplan
```

Terraform waits for the Elastic Container Service (ECS) services to reach a steady state in dev, staging and prod. The rdev root does not wait.

## Run the database migration

Happy ran the database migration after every Terraform apply. Terraform does not run it automatically. Run the migration task after the apply:

```bash
unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN

case "$TF_ROOT" in
  */prod)
    TARGET_ACCOUNT=231426846575
    CONFIG_SECRET=happy/env-prod-config
    ;;
  */stage)
    TARGET_ACCOUNT=699936264352
    CONFIG_SECRET=happy/env-stage-config
    ;;
  */dev)
    TARGET_ACCOUNT=699936264352
    CONFIG_SECRET=happy/env-dev-config
    ;;
  */rdev)
    TARGET_ACCOUNT=699936264352
    CONFIG_SECRET=happy/env-rdev-config
    ;;
  *)
    exit 1
    ;;
esac

credentials=$(AWS_PROFILE=czi-id aws sts assume-role \
  --role-arn "arn:aws:iam::${TARGET_ACCOUNT}:role/tfe-si" \
  --role-session-name single-cell-data-portal-deploy)
export AWS_ACCESS_KEY_ID=$(jq -r '.Credentials.AccessKeyId' <<<"$credentials")
export AWS_SECRET_ACCESS_KEY=$(jq -r '.Credentials.SecretAccessKey' <<<"$credentials")
export AWS_SESSION_TOKEN=$(jq -r '.Credentials.SessionToken' <<<"$credentials")

config=$(aws secretsmanager get-secret-value \
  --secret-id "$CONFIG_SECRET" \
  --query SecretString \
  --output text)
cluster=$(jq -r '.cluster_arn' <<<"$config")
subnets=$(jq -c '.private_subnets' <<<"$config")
security_groups=$(jq -c '.security_groups' <<<"$config")
task_definition=$(terraform output -raw migrate_db_task_definition_arn)

task_arn=$(aws ecs run-task \
  --cluster "$cluster" \
  --task-definition "$task_definition" \
  --launch-type FARGATE \
  --network-configuration \
    "awsvpcConfiguration={subnets=${subnets},securityGroups=${security_groups},assignPublicIp=DISABLED}" \
  --query 'tasks[0].taskArn' \
  --output text)

aws ecs wait tasks-stopped --cluster "$cluster" --tasks "$task_arn"
exit_code=$(aws ecs describe-tasks \
  --cluster "$cluster" \
  --tasks "$task_arn" \
  --query 'tasks[0].containers[0].exitCode' \
  --output text)
test "$exit_code" = "0"

unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN
unset credentials config
export AWS_PROFILE=czi-id
```

Do not continue if the migration task fails.

## Invalidate CloudFront

The frontend uses CloudFront in dev, staging and prod. Find its distribution and invalidate `index.html`:

```bash
case "$TF_ROOT" in
  */prod)
    AWS_ACCOUNT=231426846575
    DOMAIN_NAME=frontend.production.single-cell.czi.technology
    ALIAS=cellxgene.cziscience.com
    ;;
  */stage)
    AWS_ACCOUNT=699936264352
    DOMAIN_NAME=frontend.stage.single-cell.czi.technology
    ALIAS=cellxgene.staging.single-cell.czi.technology
    ;;
  */dev)
    AWS_ACCOUNT=699936264352
    DOMAIN_NAME=frontend.dev.single-cell.czi.technology
    ALIAS=cellxgene.dev.single-cell.czi.technology
    ;;
esac

credentials=$(AWS_PROFILE=czi-id aws sts assume-role \
  --role-arn "arn:aws:iam::${AWS_ACCOUNT}:role/tfe-si" \
  --role-session-name single-cell-data-portal-cloudfront)
export AWS_ACCESS_KEY_ID=$(jq -r '.Credentials.AccessKeyId' <<<"$credentials")
export AWS_SECRET_ACCESS_KEY=$(jq -r '.Credentials.SecretAccessKey' <<<"$credentials")
export AWS_SESSION_TOKEN=$(jq -r '.Credentials.SessionToken' <<<"$credentials")

distribution_id=$(aws cloudfront list-distributions \
  --query "DistributionList.Items[*].{id:Id,domain_name:Origins.Items[*].DomainName,alias:Aliases.Items[0]}[?contains(domain_name,'${DOMAIN_NAME}')&&alias=='${ALIAS}'].id" \
  --output text)
test -n "$distribution_id"
aws cloudfront create-invalidation \
  --distribution-id "$distribution_id" \
  --paths /index.html

unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN
unset credentials
export AWS_PROFILE=czi-id
```

The rdev stack does not use this CloudFront invalidation step.

## Validate

Use `terraform output` to get the frontend and backend URLs. Confirm the frontend loads, the backend health endpoint responds and the ECS services have reached their desired task counts. Run any release-specific functional or performance checks manually.

## State recovery

State lives in `s3://terragrunt-engine-state` and uses the `terragrunt-engine-state-lock` DynamoDB table. Never edit or upload the state object with the AWS CLI because that bypasses Terraform's lock-table digest.

If a state restore is necessary, stop all other Terraform work, retain a copy of the current state and use `terraform state push` with a known-good encrypted backup:

```bash
terraform state pull > current.tfstate
terraform state push backup.tfstate
```

State files contain sensitive values. Store backups only in an approved encrypted location and delete local copies when recovery is complete.
