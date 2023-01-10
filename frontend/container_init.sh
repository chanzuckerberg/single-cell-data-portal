#!/bin/bash
echo
echo " ====="
echo "| starting frontend container"
echo " ====="
echo

# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

# no verbose
set -x # config

# We build the next.js app during our CI process, but then we need to be able
# to run the resulting docker image in dev/rdev/staging/prod. This isn't really
# possible to control well with next.js native utilities, so we have to rewrite
# any environment-specific configuration in our *pre-built* nextjs app on
# container startup.
#
# This function is *intentionally* idempotent - it can be executed multiple
# times with different api_url / deployment_env vars and it should get the
# built app to a sane state on each run.
function apply_path {
  envDataOriginalFile="./runtime_configs/${DEPLOYMENT_STAGE}.env"
  if [ ! -e $envDataOriginalFile ]; then
    echo "ERROR: Invalid DEPLOYMENT_STAGE var, please fix me"
    exit 1
  fi
  envDataFile="./runtime_configs/current_config"
  nextFolder='.next'
  nextFolderBackup='.next.original'

  # Expand any env vars in the original runtime var file and write it
  # to the file we're going to use for replacing vars in our next build.
  envsubst < $envDataOriginalFile > $envDataFile

  # Make a backup of the original next build files on the first run
  # and in subsequent runs, start from that backup set of files
  # so we can safely run this script multiple times.
  if [ ! -e $nextFolderBackup ]; then
    cp -r $nextFolder $nextFolderBackup
  else
    cp -r $nextFolderBackup/. $nextFolder
  fi

  # read in the runtime config file a line at a time, reading variable names and values.
  # The variable names need to match the names we've given in our src/configs/build.js file!
  # Then any occurrences of these var names are replaced with their values in the .next dir.
  while read line; do
    # no comment or not empty
    if [ "${line:0:1}" == "#" ] || [ "${line}" == "" ]; then
      continue
    fi

    # read variable names and values
    configName="$(cut -d'=' -f1 <<<"$line")"
    configValue="$(cut -d'=' -f2 <<<"$line")"    # get system env

    # replace the variable names in our built app with the appropriate values. Use rg (ripgrep) for speed.
    echo "Replace: ${configName} with: ${configValue}"
    rg $configName --files-with-matches -0 $nextFolder | xargs -0 sed -i "s#$configName#$configValue#g"
  done < $envDataFile
}

# Build and run without dev mode in remote dev env.
if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  cp src/configs/local.js src/configs/configs.js
  mkdir -p ./node_modules/.next-dev-mobile
  cp -r /tmp/pkcs12/* ./node_modules/.next-dev-mobile
  # next-dev-https looks for cert.pem and key.pem in node_modules/.next-dev-mobile/ (otherwise generates anew 👎)
  # We want the dev server to find this key and cert because they've already been added to our keychain
  mv ./node_modules/.next-dev-mobile/server.crt ./node_modules/.next-dev-mobile/cert.pem
  mv ./node_modules/.next-dev-mobile/server.key ./node_modules/.next-dev-mobile/key.pem
  exec npm run dev
else
  # We need "-- --" because `npm run build-and-start-prod`
  # runs `npm run build && npm run serve` under the hood,
  # so we need to pass `-- -p 9000` to `npm run serve`, which
  # will then call `next start -p 9000` correctly
  apply_path
  exec npm run serve -- -- -p 9000
fi
