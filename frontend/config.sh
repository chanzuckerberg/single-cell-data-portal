#!/bin/bash
if [ ! -z "$API_URL" ]; then
  cat src/configs/$DEPLOYMENT_STAGE.js | envsubst >src/configs/configs.js
else
  cp src/configs/local.js src/configs/configs.js
fi
