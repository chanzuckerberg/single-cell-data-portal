#!/bin/bash

checks=$(curl -s  "https://api.github.com/repos/chanzuckerberg/single-cell-data-portal/actions/runs?actor=czibuildbot" | jq -r '.workflow_runs[] | select(.display_title == "stage") | .conclusion' | head -n 1)
if [ "$checks" == "success" ]; then
  echo "All checks passed in most recent staging deployment!"
else
  echo "Some checks failed in most recent staging deployment. Please check GHA runs."
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
git reset --hard origin staging
echo "Checking out 'prod' branch and pull latest"
git checkout prod
git reset --hard origin prod
echo "Confirming checked out branch receiving merge is: $(git branch --show-current)"

echo "Latest commit on 'prod' branch is: $(git log -1 --format='%H' prod)"
echo "Latest commit on 'staging' branch is: $(git log -1 --format='%H' staging)"

echo "About to merge 'staging' branch into 'prod' branch"
git merge --verbose staging -m "Merging staging branch into prod branch"
echo "Merge completed"

echo "Pushing to Prod"
git push origin prod