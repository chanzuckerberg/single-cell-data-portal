#!/usr/bin/env bash

set -euo pipefail

service="${1:-}"
repository="${2:-}"
profile="${3:-}"
image_tag="${4:-}"
branch_tag="${5:-}"
commit_sha="${6:-}"
branch_name="${7:-}"

[[ -n "$service" && -n "$repository" && -n "$profile" && -n "$image_tag" && -n "$branch_tag" && -n "$commit_sha" && -n "$branch_name" ]] || {
  echo "Usage: scripts/build-images.sh <service> <repository> <profile> <image-tag> <branch-tag> <commit-sha> <branch-name>" >&2
  exit 1
}
[[ -n "${DOCKER_REPO:-}" ]] || {
  echo "DOCKER_REPO is required" >&2
  exit 1
}

export HAPPY_COMMIT="$commit_sha"
export HAPPY_BRANCH="$branch_name"
export HAPPY_TAG="$image_tag"

docker compose --profile "$profile" build "$service"

source_image="${DOCKER_REPO}${repository}:latest"
for tag in "$image_tag" "$branch_tag"; do
  target_image="${DOCKER_REPO}${repository}:${tag}"
  docker tag "$source_image" "$target_image"
  docker push "$target_image"
done
