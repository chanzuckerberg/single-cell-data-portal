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
cp src/configs/local.js src/configs/configs.js
exec gatsby develop --host 0.0.0.0
