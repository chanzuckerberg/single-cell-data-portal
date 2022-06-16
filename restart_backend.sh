#!/bin/bash
docker restart $(docker ps -f name=backend -q)
./scripts/setup_dev_data.sh