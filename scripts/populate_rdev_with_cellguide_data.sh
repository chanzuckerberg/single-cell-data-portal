#!/bin/bash
STACK_NAME=""
SRC_DEPLOYMENT="staging"
POPULATE_ONLY_GPT=false

function usage {
    echo "Usage: $0 [-s|--src_deployment <src_deployment>] [-g|--gpt_only] <stack_name>"
    echo "Make sure AWS_PROFILE is set to single-cell-dev prior to running this command."
    echo ""
    echo "Populates the specified rdev stack with data."
    echo ""
    echo "Options:"
    echo "  -s, --src_deployment    The source deployment to use. Can be either 'staging' or 'dev'. Default is 'staging'."
    echo "  -g, --gpt_only          Only populate gpt_descriptions and gpt_seo_descriptions. Default is false."
    echo ""
    echo "Arguments:"
    echo "  stack_name              The name of the rdev stack to be populated."
    exit 1
}

while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage
      ;;
    -s|--src_deployment)
      SRC_DEPLOYMENT="$2"
      if [[ "$SRC_DEPLOYMENT" != "staging" && "$SRC_DEPLOYMENT" != "dev" ]]; then
        echo "Error: SRC_DEPLOYMENT must be either 'staging' or 'dev'"
        exit 1
      fi
      shift 2
      ;;
    -g|--gpt_only)
      POPULATE_ONLY_GPT=true
      shift
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      STACK_NAME="$1"
      shift
      ;;
  esac
done

if [ -z "$STACK_NAME" ]
then
  echo "Error: Stack name is required"
  exit 1
fi

LATEST_SNAPSHOT_IDENTIFIER=$(aws s3 cp s3://cellguide-data-public-${SRC_DEPLOYMENT}/latest_snapshot_identifier -)
if [ -n "$LATEST_SNAPSHOT_IDENTIFIER" ]
then
  aws s3 sync s3://cellguide-data-public-${SRC_DEPLOYMENT}/gpt_seo_descriptions s3://cellguide-data-public-dev/env-rdev-cellguide/${STACK_NAME}/gpt_seo_descriptions
  aws s3 sync s3://cellguide-data-public-${SRC_DEPLOYMENT}/gpt_descriptions s3://cellguide-data-public-dev/env-rdev-cellguide/${STACK_NAME}/gpt_descriptions
  if [ "$POPULATE_ONLY_GPT" = false ] ; then
    aws s3 sync s3://cellguide-data-public-${SRC_DEPLOYMENT}/${LATEST_SNAPSHOT_IDENTIFIER} s3://cellguide-data-public-dev/env-rdev-cellguide/${STACK_NAME}/${LATEST_SNAPSHOT_IDENTIFIER}
    aws s3 cp s3://cellguide-data-public-${SRC_DEPLOYMENT}/latest_snapshot_identifier s3://cellguide-data-public-dev/env-rdev-cellguide/${STACK_NAME}/latest_snapshot_identifier
  fi
else
  echo "Warning: LATEST_SNAPSHOT_IDENTIFIER is invalid"
fi

