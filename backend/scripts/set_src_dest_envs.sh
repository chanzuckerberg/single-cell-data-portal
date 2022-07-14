#!/usr/bin/env bash

export SRC_ENV=${1:-prod}
export DEST_ENV=${2:-dev}

if [[ ! "$SRC_ENV" =~ ^(dev|staging|prod)$ ]]; then
  echo ERROR: invalid value for source deployment env: \"$SRC_ENV\"
  exit 1
fi

if [[ ! "$DEST_ENV" =~ ^(dev|staging)$ ]]; then
  echo ERROR: invalid value for destination deployment env: \"$DEST_ENV\"
  exit 1
fi

if [[ "$DEST_ENV" == staging && "$SRC_ENV" == dev ]]; then
  echo ERROR: cannot mirror from $SRC_ENV to $DEST_ENV
  exit 1
fi
