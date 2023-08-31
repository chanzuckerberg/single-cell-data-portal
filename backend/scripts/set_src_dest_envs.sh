#!/usr/bin/env bash

export SRC_ENV=${1:-prod}
export DEST_ENV=${2:-dev}

# require DEST_ENV=rdev when STACK is set
if [[ ! -z "$STACK" ]]; then
    if [[ ! "$DEST_ENV" == rdev ]]; then
        echo "ERROR: when STACK var is set, DEST_ENV must be 'rdev'"
        exit 1
    fi
fi

# allow dev->dev for testing
if [[ "$SRC_ENV" == dev && "$DEST_ENV" == dev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the dev env to itself. Useful for testing the mirroring operation, but otherwise useless. Proceed? (Y/n) " ANS
        [[ $ANS == 'Y' ]] || exit 1
    fi
# allow dev->rdev for testing
elif [[ "$SRC_ENV" == dev && "$DEST_ENV" == rdev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the dev env to rdev. Useful for testing fixes in rdev for data-related bugs found in dev, but otherwise useless. Proceed? (Y/n) " ANS
        [[ $ANS == 'Y' ]] || exit 1
    fi
# allow staging->dev for testing
elif [[ "$SRC_ENV" == staging && "$DEST_ENV" =~ (dev|rdev) ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the staging env to ${DEST_ENV}. Useful for testing fixes in ${DEST_ENV} for data-related bugs found in staging, but otherwise useless. Proceed? (Y/n) " ANS
        [[ $ANS == 'Y' ]] || exit 1
    fi
# allow prod->dev and prod->staging
elif [[ ! ("$SRC_ENV" == prod && "$DEST_ENV" =~ ^(dev|rdev|staging)$) ]]; then
    echo "ERROR: invalid value src/dest envs: $SRC_ENV -> $DEST_ENV"
    exit 1
fi


