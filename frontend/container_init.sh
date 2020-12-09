#!/bin/bash
echo
echo " ====="
echo "| starting frontend container"
echo " ====="
echo
if [ ! -e ./node_modules ]; then
  mkdir -p ./node_modules
fi
ln -sf /opt/node_app/node_modules/* ./node_modules/.
if [ ! -z "$API_URL" ]; then
  cat src/configs/rdev.js | envsubst > src/configs/configs.js
else
  cp src/configs/local.js src/configs/configs.js
fi
exec gatsby develop --host 0.0.0.0
