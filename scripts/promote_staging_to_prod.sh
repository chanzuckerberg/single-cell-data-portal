#!/bin/bash

if [ ! -z $1 ]; then
  most_recent_deploy_run=$(curl -s -H "Authorization: Bearer $1" "https://api.github.com/repos/chanzuckerberg/single-cell-data-portal/actions/workflows/deploy-happy-stack.yml/runs" | jq -r '[.workflow_runs[] | select(.display_title == "stage")][0]')
else
  most_recent_deploy_run=$(curl -s  "https://api.github.com/repos/chanzuckerberg/single-cell-data-portal/actions/workflows/deploy-happy-stack.yml/runs" | jq -r '[.workflow_runs[] | select(.display_title == "stage")][0]')
fi

checks=$(echo "$most_recent_deploy_run" | jq -r '.conclusion')
gha_head_sha=$(echo "$most_recent_deploy_run" | jq -r '.head_sha')

if [ "$checks" == "success" ]; then
  echo "All checks passed in most recent staging deployment!"
else
  echo "Not all checks have concluded successfully in most recent staging deployment. Please check GHA runs."
  exit 1
fi

if [[ $(git diff HEAD) ]]; then
  echo "Local has uncommitted changes, please commit or stash and try again."
  exit 1
fi

echo "fetch origin branch history"
git fetch origin
echo "Checking out 'staging' branch and pull latest"
git checkout staging
git reset --hard origin/staging
echo "Checking out 'prod' branch and pull latest"
git checkout prod
git reset --hard origin/prod
echo "Confirming checked out branch receiving merge is: $(git branch --show-current)"

# Get most recent commit, excluding commits by 'GitHub Actions' (author for merge commits that do not trigger deployments)
prod_head_sha=$(git log -1 --format='%H' prod --perl-regexp --author='^(?!(.*(GitHub Actions)))')
staging_head_sha=$(git log -1 --format='%H' staging --perl-regexp --author='^(?!(.*(GitHub Actions)))')

echo "Latest commit on 'prod' branch is: $prod_head_sha"
echo "Latest commit on 'staging' branch is: $staging_head_sha"

if [ "$staging_head_sha" != "$gha_head_sha" ]; then
  echo "Latest commit on 'staging' branch ($staging_head_sha) does not match commit for staging deployment GHA checks ($gha_head_sha).
  A new commit may have been introduced to staging, please check repo and try running again."
  exit 1
fi

echo "About to merge 'staging' branch into 'prod' branch"
if git merge --verbose staging -m "Merging staging branch into prod branch"; then
  echo "Merge was Successful"
else
  echo "Merge has conflicts or other issues. Please resolve and try again."
  exit 1
fi

echo "Pushing to Prod"
if git push origin prod; then
  echo "Successfully Promoted Staging to Prod"
else
  echo "Staging push to Prod failed"
  exit 1
fi