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

# Build and run without dev mode in remote dev env.
if [ "${DEPLOYMENT_STAGE}" == "test" ]; then
  exec npm run dev
else
  exec npm run serve -p 9000
fi
