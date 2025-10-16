#!/usr/bin/env bash

echo "Setting env vars for src/dest envs"
export SRC_ENV=${1:-prod}
export DEST_ENV=${2:-dev}
export SRC_PORT=${3:-5432}
export DEST_PORT=${4:-5433}

# DRY_RUN can be passed as the 10th parameter or as an environment variable
if [[ -n "${10}" ]] || [[ -n "$DRY_RUN" ]]; then
    export DRY_RUN=1
    echo "🔍 DRY-RUN MODE: No changes will be made"
fi

# Forbid rdev as SRC_ENV
if [[ "$SRC_ENV" == rdev ]]; then
  echo "ERROR: mirroring FROM rdev is not supported; only mirroring TO rdev is supported."
  exit 1
fi

# require DEST_ENV=rdev when STACK is set
if [[ -n "$STACK" ]]; then
    if [[ "$DEST_ENV" != rdev ]]; then
        echo "ERROR: when STACK var is set, DEST_ENV must be 'rdev'"
        exit 1
    fi
fi

# allow dev->dev for testing
if [[ "$SRC_ENV" == dev && "$DEST_ENV" == dev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the dev env to itself. Useful for testing the mirroring operation, but otherwise useless. Proceed? (Y/n) " ANS
        echo
        [[ $ANS == 'Y' ]] || exit 1
        export NO_PROMPT=1  # Don't prompt again in subsequent scripts
    fi
# allow dev->rdev for testing
elif [[ "$SRC_ENV" == dev && "$DEST_ENV" == rdev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the dev env to rdev. Useful for testing fixes in rdev for data-related bugs found in dev, but otherwise useless. Proceed? (Y/n) " ANS
        echo
        [[ $ANS == 'Y' ]] || exit 1
        export NO_PROMPT=1  # Don't prompt again in subsequent scripts
    fi
# allow staging->dev for testing
elif [[ "$SRC_ENV" == staging && "$DEST_ENV" =~ (dev|rdev) ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the staging env to ${DEST_ENV}. Useful for testing fixes in ${DEST_ENV} for data-related bugs found in staging, but otherwise useless. Proceed? (Y/n) " ANS
        echo
        [[ $ANS == 'Y' ]] || exit 1
        export NO_PROMPT=1  # Don't prompt again in subsequent scripts
    fi
# allow prod->dev and prod->staging
elif [[ ! ("$SRC_ENV" == prod && "$DEST_ENV" =~ ^(dev|rdev|staging)$) ]]; then
    echo "ERROR: invalid value src/dest envs: $SRC_ENV -> $DEST_ENV"
    exit 1
fi


