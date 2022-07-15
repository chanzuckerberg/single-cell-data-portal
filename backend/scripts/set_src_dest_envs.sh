#!/usr/bin/env bash

export SRC_ENV=${1:-prod}
export DEST_ENV=${2:-dev}

# allow prod->dev and prod->staging
if [[ "$SRC_ENV" == prod && "$DEST_ENV" =~ ^(dev|staging)$ ]]; then
    exit 0
fi

# allow dev->dev for testing
if [[ "$SRC_ENV" == dev && "$DEST_ENV" == dev ]]; then
    read -n 1 -p "ATTENTION: You are mirroring the dev env to itself. Useful for testing the mirroring operation, but otherwise useless. Proceed? (Y/n) " ANS
    [[ $ANS == 'Y' ]] || exit 1
    exit 0
fi

# allow staging->dev for testing
if [[ "$SRC_ENV" == staging && "$DEST_ENV" == dev ]]; then
    read -n 1 -p "ATTENTION: You are mirroring the staging env to dev. Useful for testing fixes in dev for data-related bugs found in staging, but otherwise useless. Proceed? (Y/n) " ANS
    [[ $ANS == 'Y' ]] || exit 1
    exit 0
fi

echo "ERROR: invalid value src/dest envs: $SRC_ENV -> $DEST_ENV"
exit 1
