# Manual deployment

The data portal no longer deploys through Happy or Terraform Enterprise (TFE). A core infrastructure engineer must use `scripts/deploy.sh` for a break-glass deployment.

## Supported deployments

The script supports three deployment paths:

| Target  | Source                    | Terraform root                | Stack        |
| ------- | ------------------------- | ----------------------------- | ------------ |
| rdev    | An open pull request      | `.happy/terraform/envs/rdev`  | `rdevstack`  |
| staging | The current `main` commit | `.happy/terraform/envs/stage` | `stagestack` |
| prod    | The current `main` commit | `.happy/terraform/envs/prod`  | `prodstack`  |

The rdev root manages one shared stack. It does not create a stack for each pull request, and closing a pull request does not remove it.

## Prerequisites

- Terraform 1.3
- authenticated `aws`, `gh`, `git` and `jq` command-line tools
- AWS access through the `czi-id` profile
- permission to assume `tfe-si` in the target account and `core-platform-prod`
- a clean checkout of the exact commit to deploy

Run the script from the repository root. It stops if the working tree is dirty, the checkout does not match the requested source or a fixed environment has empty state.

## Deploy a pull request to rdev

Use rdev to test an open pull request before merge:

```bash
gh pr checkout 1234
scripts/deploy.sh rdev --pr 1234
```

The script confirms that the local commit matches the pull request head. Cross-repository pull requests are not supported because GitHub cannot dispatch the image workflow against a fork branch.

Coordinate with other operators before using rdev. A deployment replaces the application version in the shared `rdevstack`.

## Deploy main to staging

Update the local `main` branch to exactly match `origin/main`, then deploy:

```bash
git checkout main
git pull --ff-only origin main
scripts/deploy.sh staging
```

## Deploy main to prod

Production also deploys from `main`. Update the local branch, review the commit and run:

```bash
git checkout main
git pull --ff-only origin main
scripts/deploy.sh prod
```

The script prompts before applying the reviewed Terraform plan.

## What the script does

For every target, `scripts/deploy.sh`:

1. Resolves the exact pull request or `main` commit and verifies the local checkout.
2. Checks all nine ECR repositories for the corresponding `sha-<first eight characters>` image tag.
3. Reuses the images when every repository already has the tag. Otherwise, it triggers the `Build Images` GitHub Actions workflow and waits for every image build job to succeed.
4. Verifies that the workflow built the requested commit.
5. Initializes the target Terraform root and creates a saved plan using that image tag.
6. Displays the plan and asks for confirmation.
7. Applies the saved plan.
8. Runs the database migration task and verifies its exit code.
9. Invalidates CloudFront for staging and prod.
10. Prints the Terraform outputs for validation.

The image workflow uses Docker Compose to build every environment's images and pushes them to the development Elastic Container Registry (ECR) repositories. Staging and prod also pull from those repositories. The nine images build in parallel and import inline BuildKit cache from their `branch-main` images. Image building does not use Happy, TFE or Terraform.

Do not derive the image tag from a local commit. The script uses the workflow run's actual head SHA and aborts if it differs from the requested deployment commit.

Cancelling at the Terraform apply prompt does not discard built images. Run the same deployment command again from the same commit. The script finds all nine SHA-tagged images in ECR, skips the build and returns to Terraform planning.

If the branch advanced after the images were built, pass the previous tag explicitly:

```bash
scripts/deploy.sh rdev --pr 1234 --image-tag sha-01234567
```

The script verifies that all nine images exist before using an explicit tag. Use this option only when you intend to deploy images from a commit other than the current checkout.

## Review the plan

Stop if a fixed environment plans to create every resource. That indicates an empty or incorrect state key.

Review every replacement and deletion before answering the apply prompt. The dev state predates the current configuration and contains substantial drift, but this script does not support deploying dev. Staging and prod can also contain drift that is unrelated to the application commit.

## Validate the deployment

Use the outputs printed by the script to find the frontend and backend URLs. Confirm that:

- the frontend loads
- the backend health endpoint responds
- the Elastic Container Service (ECS) services reach their desired task counts
- release-specific functional or performance checks pass

## State recovery

State lives in `s3://terragrunt-engine-state` and uses the `terragrunt-engine-state-lock` DynamoDB table. Never edit or upload a state object with the AWS CLI because that bypasses Terraform's lock-table digest.

If a state restore is necessary, stop all other Terraform work, retain a copy of the current state and use `terraform state push` with a known-good encrypted backup:

```bash
terraform state pull > current.tfstate
terraform state push backup.tfstate
```

State files contain sensitive values. Store backups only in an approved encrypted location and delete local copies when recovery is complete.
