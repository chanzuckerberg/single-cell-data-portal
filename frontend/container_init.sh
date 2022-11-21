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
  cat src/configs/$DEPLOYMENT_STAGE.js | envsubst >src/configs/configs.js
else
  cp src/configs/local.js src/configs/configs.js
fi

# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

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
  exec npm run build-and-start-prod -- -- -p 9000
fi
