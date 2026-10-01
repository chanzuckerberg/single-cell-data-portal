#!/usr/bin/env bash

set -euo pipefail

slice="${1:-}"
image_tag="${2:-}"
branch_tag="${3:-}"
commit_sha="${4:-}"
branch_name="${5:-}"

[[ -n "$slice" && -n "$image_tag" && -n "$branch_tag" && -n "$commit_sha" && -n "$branch_name" ]] || {
  echo "Usage: scripts/build-images.sh <slice> <image-tag> <branch-tag> <commit-sha> <branch-name>" >&2
  exit 1
}
[[ -n "${DOCKER_REPO:-}" ]] || {
  echo "DOCKER_REPO is required" >&2
  exit 1
}

case "$slice" in
  frontend)
    profile=fullstack
    services=(frontend backend backend-de backend-wmg)
    repositories=(corpora-frontend corpora-backend corpora-backend-de corpora-backend-wmg)
    ;;
  upload_failures)
    profile=upload_failures
    services=(upload_failures)
    repositories=(corpora-upload-failures)
    ;;
  upload_success)
    profile=upload_success
    services=(upload_success)
    repositories=(corpora-upload-success)
    ;;
  processing)
    profile=processing
    services=(processing)
    repositories=(corpora-upload)
    ;;
  wmg_processing)
    profile=wmg_processing
    services=(wmg_processing)
    repositories=(wmg-processing)
    ;;
  cellguide_pipeline)
    profile=cellguide_pipeline
    services=(cellguide_pipeline)
    repositories=(cellguide-pipeline)
    ;;
  *)
    echo "Unsupported image slice: $slice" >&2
    exit 1
    ;;
esac

export HAPPY_COMMIT="$commit_sha"
export HAPPY_BRANCH="$branch_name"
export HAPPY_TAG="$image_tag"

for repository in "${repositories[@]}"; do
  docker pull "${DOCKER_REPO}${repository}:branch-main" || true
done

docker compose --profile "$profile" build "${services[@]}"

for repository in "${repositories[@]}"; do
  source_image="${DOCKER_REPO}${repository}:latest"
  for tag in "$image_tag" "$branch_tag"; do
    target_image="${DOCKER_REPO}${repository}:${tag}"
    docker tag "$source_image" "$target_image"
    docker push "$target_image"
  done
done
