#!/bin/bash
echo
echo " ====="
echo "| starting frontend container"
echo " ====="
echo
if [ ! -e ./node_modules ]; then
  mkdir -p ./node_modules
fi
ln -sf /opt/node_app/node_modules/* /opt/node_app/node_modules/.bin ./node_modules/.
if [ ! -z "$API_URL" ]; then
  envsubst < runtime_configs/$DEPLOYMENT_STAGE.env | sponge src/configs/configs.js
else
  cp src/configs/local.js src/configs/configs.js
fi


# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

# no verbose
set -x # config

function apply_path {
  envDataFile="./runtime_configs/${DEPLOYMENT_STAGE}.env"
  nextFolder='.next'
  nextFolderBackup='.next.original'

  # Make a backup of the original next build files on the first run
  # and in subsequent runs, start from that original set of files
  # so we can safely run this script multiple times.
  if [ ! -e $nextFolderBackup ]; then
    cp -r $nextFolder $nextFolderBackup
  else
    cp -r $nextFolderBackup/. $nextFolder
  fi

  # read all config file  
  while read line; do
    # no comment or not empty
    if [ "${line:0:1}" == "#" ] || [ "${line}" == "" ]; then
      continue
    fi
    
    # split
    configName="$(cut -d'=' -f1 <<<"$line")"
    configValue="$(cut -d'=' -f2 <<<"$line")"    # get system env
    
    # replace all
    echo "Replace: ${configName} with: ${configValue}"
    rg $configName --files-with-matches -0 $nextFolder | xargs -0 sed -i "s#$configName#$configValue#g"
  done < $envDataFile
}

# Build and run without dev mode in remote dev env.
if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
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
