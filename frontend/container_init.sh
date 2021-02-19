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
  cat src/configs/rdev.js | envsubst > src/configs/configs.js
else
  cp src/configs/local.js src/configs/configs.js
fi

# If user passed a command line, run it in place of the server
if [ $# -ne 0 ]; then
  exec "$@"
fi

# Build and run without dev mode in remote dev env.
if [ "${DEPLOYMENT_STAGE}" == "dev" ]; then
    exec gatsby develop --host 0.0.0.0
else
    npm run build
    exec gatsby serve --host 0.0.0.0
fi
