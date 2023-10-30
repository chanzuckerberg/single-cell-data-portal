#!/bin/bash
STACK_NAME=""
SRC_DEPLOYMENT="staging"
SNAPSHOT_VERSION="v3"

function usage {
  echo "Usage: $0 [-s|--src_deployment <src_deployment>] [-g|--gpt_only] <stack_name>"
  echo "Make sure AWS_PROFILE is set to single-cell-dev prior to running this command."
  echo ""
  echo "Populates the specified rdev stack with data."
  echo ""
  echo "Options:"
  echo "  -s, --src_deployment    The source deployment to use. Can be either 'staging' or 'dev'. Default is 'staging'."
  echo "  -v, --version           The snapshot schema version to mirror. Default is 'v1'."
  echo ""
  echo "Arguments:"
  echo "  stack_name              The name of the rdev stack to be populated."
  exit 1
}

while (("$#")); do
  case "$1" in
  -h | --help)
    usage
    ;;
  -s | --src_deployment)
    SRC_DEPLOYMENT="$2"
    if [[ "$SRC_DEPLOYMENT" != "staging" && "$SRC_DEPLOYMENT" != "dev" ]]; then
      echo "Error: SRC_DEPLOYMENT must be either 'staging' or 'dev'"
      exit 1
    fi
    shift 2
    ;;
  -v | --version)
    SNAPSHOT_VERSION="$2"
    shift 2
    ;;
  --) # end argument parsing
    shift
    break
    ;;
  -* | --*=) # unsupported flags
    echo "Error: Unsupported flag $1" >&2
    exit 1
    ;;
  *) # preserve positional arguments
    STACK_NAME="$1"
    shift
    ;;
  esac
done

if [ -z "$STACK_NAME" ]; then
  echo "Error: Stack name is required"
  exit 1
fi

LATEST_SNAPSHOT_IDENTIFIER=$(aws s3 cp s3://cellxgene-wmg-${SRC_DEPLOYMENT}/snapshots/${SNAPSHOT_VERSION}/latest_snapshot_identifier -)
if [ -n "$LATEST_SNAPSHOT_IDENTIFIER" ]; then
  aws s3 rm s3://env-rdev-wmg/${STACK_NAME}/snapshots/${SNAPSHOT_VERSION}/${LATEST_SNAPSHOT_IDENTIFIER} --recursive
  aws s3 sync s3://cellxgene-wmg-${SRC_DEPLOYMENT}/snapshots/${SNAPSHOT_VERSION}/${LATEST_SNAPSHOT_IDENTIFIER} s3://env-rdev-wmg/${STACK_NAME}/snapshots/${SNAPSHOT_VERSION}/${LATEST_SNAPSHOT_IDENTIFIER}

  # Check the exit code of the last command (AWS CLI)
  if [ $? -ne 0 ]; then
    echo "Error: AWS CLI command failed"
    exit 1
  fi

  aws s3 cp s3://cellxgene-wmg-${SRC_DEPLOYMENT}/snapshots/${SNAPSHOT_VERSION}/latest_snapshot_identifier s3://env-rdev-wmg/${STACK_NAME}/snapshots/${SNAPSHOT_VERSION}/latest_snapshot_identifier

  # Check the exit code of the last command (AWS CLI)
  if [ $? -ne 0 ]; then
    echo "Error: AWS CLI command failed"
    exit 1
  fi
else
  echo "Warning: LATEST_SNAPSHOT_IDENTIFIER is invalid"
  exit 1
fi
