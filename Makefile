SHELL:=/bin/bash

.PHONY: fmt
fmt:
	black .

.PHONY: lint
lint:
	flake8 backend tests

.PHONY: unit-test
unit-test:
	-docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres
	PYTHONWARNINGS=ignore:ResourceWarning python3 -m coverage run \
		-m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose; \
	test_result=$$?; \
	$(MAKE) clean_test_db; \
	exit $$test_result


clean_test_db:
	-docker stop test_db
	-docker rm test_db

.PHONY: local-server
local-server: local-backend local-frontend

.PHONY: functional-test
functional-test:
	python3 -m unittest discover --start-directory tests/functional --top-level-directory . --verbose

.PHONY: local-database
local-database: clean_test_db
	docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres
	python ./scripts/populate_db.py


.PHONY: local-backend
local-backend: local-database
	$(MAKE) local-server -C ./backend/chalice/api_server DEPLOYMENT_STAGE=test

.PHONY: local-frontend
local-frontend:
	$(MAKE) local-server -C ./frontend DEPLOYMENT_STAGE=test
