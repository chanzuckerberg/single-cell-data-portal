#!/usr/bin/env bash

cd ../
docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres
make package -C ./backend/chalice/api_server
python run_local_server.py &
cd ../frontend
npm install
gatsby develop
