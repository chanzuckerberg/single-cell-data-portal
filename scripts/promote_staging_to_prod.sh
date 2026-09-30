#!/bin/bash

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