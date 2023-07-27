#!/usr/bin/env bash

export SRC_ENV=${1:-prod}
export DEST_ENV=${2:-dev}
export RDEV_ENV=${3:-}

# allow dev->dev for testing
if [[ "$SRC_ENV" == dev && "$DEST_ENV" == dev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the dev env to itself. Useful for testing the mirroring operation, but otherwise useless. Proceed? (Y/n) " ANS
        [[ $ANS == 'Y' ]] || exit 1
    fi
# allow staging->dev for testing
elif [[ "$SRC_ENV" == staging && "$DEST_ENV" == dev ]]; then
    if [[ -z $NO_PROMPT ]]; then
        read -n 1 -p "ATTENTION: You are mirroring the staging env to dev. Useful for testing fixes in dev for data-related bugs found in staging, but otherwise useless. Proceed? (Y/n) " ANS
        [[ $ANS == 'Y' ]] || exit 1
    fi
# allow dev||staging->rdev for creating test environments
elif [[ ("$SRC_ENV" =~ ^(dev|staging)$ && "$DEST_ENV" == rdev) ]]; then
  if [[ -z "$RDEV_ENV" ]]; then
    echo "ERROR: You must specify an rdev environment name when mirroring dev->rdev"
    exit 1
  elif [[ -z $NO_PROMPT ]]; then
      read -n 1 -p "ATTENTION: You are mirroring the ${SRC_ENV} env to rdev/${RDEV_ENV}. Useful for testing fixes in dev for
      data-related bugs found in staging, but otherwise useless. Proceed? (Y/n) " ANS
      [[ $ANS == 'Y' ]] || exit 1
  fi
# allow prod->dev and prod->staging
elif [[ ! ("$SRC_ENV" == prod && "$DEST_ENV" =~ ^(dev|staging)$) ]]; then
    echo "ERROR: invalid value src/dest envs: $SRC_ENV -> $DEST_ENV"
    exit 1
fi


