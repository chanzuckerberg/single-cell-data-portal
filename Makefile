SHELL:=/bin/bash
  # Force using BuildKit instead of normal Docker, required so that metadata
  # is written/read to allow us to use layers of previous builds as cache.
export DOCKER_BUILDKIT:=1
export COMPOSE_DOCKER_CLI_BUILD:=1
export COMPOSE_OPTS:=--env-file .env.ecr
ifeq ($(AWS_ACCESS_KEY_ID),)
	export TEST_AWS_PROFILE ?= single-cell-dev
endif
export DEPLOYMENT_STAGE ?= test
COVERAGE_DATA_FILE=.coverage.$(shell git rev-parse --short HEAD)
COVERAGE_FILTER=--omit=*/dist-packages/*,*/backend/database/*,*/backend/scripts/*,*/site-packages/* --include=*/backend/*,*/tests/unit/*
export COVERAGE_RUN_ARGS:=--data-file=$(COVERAGE_DATA_FILE) --parallel-mode $(COVERAGE_FILTER)

.PHONY: fmt
fmt:
	black --config=pyproject.toml backend scripts tests
	ruff check --fix --config=pyproject.toml backend tests scripts

.PHONY: lint
lint:
	ruff check --config=pyproject.toml backend tests

.PHONY: unit-test
unit-test: local-unit-test
	# Keeping old target name for reverse comatibility

.PHONY: wmg-processing-unittest
wmg-processing-unittest:
	# This target is intended to be run INSIDE the wmg processing container
	DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning coverage run $(COVERAGE_RUN_ARGS) -m pytest \
	tests/unit/wmg_processing/ --rootdir=. --alluredir=./allure-results --verbose;

.PHONY: cellguide-pipeline-unittest
cellguide-pipeline-unittest:
	# This target is intended to be run INSIDE the cellguide pipeline container
	DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning coverage run $(COVERAGE_RUN_ARGS) -m pytest \
	tests/unit/cellguide_pipeline/ --rootdir=. --alluredir=./allure-results --verbose;

.PHONY: functional-test
functional-test: local-functional-test
	# Keeping old target name for reverse compatibility

.PHONY: container-functionaltest
container-functionaltest:
	# This target is intended to be run INSIDE a container
	python3 -m unittest discover --start-directory tests/functional/ --top-level-directory . --verbose

.PHONY: prod-performance-test
prod-performance-test:
	python3 -m unittest discover --start-directory tests/performance --top-level-directory . --verbose

.PHONY: e2e
e2e:
	$(MAKE) e2e -C ./frontend DEPLOYMENT_STAGE=$(DEPLOYMENT_STAGE)

help: ## display help for this makefile
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help

# ------ LOCAL-DEV ------
# Local-dev related commands are below this line for now.
oauth/pkcs12/certificate.pfx:
	# All calls to the openssl cli happen in the oidc-server-mock container.
	@echo "Generating certificates for local dev"
	docker run -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint bash soluto/oidc-server-mock:0.3.0 ./generate_cert.sh
	@if [ "$$(uname -s)" == "Darwin" ]; then \
		echo "Installing generated certs into the local keychain requires sudo access:"; \
		sudo security add-trusted-cert -d -p ssl -k /Library/Keychains/System.keychain oauth/pkcs12/server.crt; \
	fi
	# Linux assumes Ubuntu
	if [ "$$(uname -s)" == "Linux" ]; then \
		sudo cp oauth/pkcs12/server.crt /usr/local/share/ca-certificates/; \
		sudo update-ca-certificates; \
	fi
	docker run -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint bash soluto/oidc-server-mock:0.3.0 ./generate_pfx.sh
	# On Linux, the pkcs12 directory gets written to with root permission. Force ownership to our user.
	sudo chown -R $$(id -u):$$(id -g) $(PWD)/oauth/pkcs12

.PHONY: .env.ecr
.env.ecr:
	echo DOCKER_REPO=$$(aws sts get-caller-identity --profile single-cell-dev | jq -r .Account).dkr.ecr.us-west-2.amazonaws.com/ > .env.ecr;
	echo "HAPPY_COMMIT=$(shell git rev-parse --verify HEAD)" >> .env.ecr
	echo "HAPPY_BRANCH=$(shell git branch --show-current)" >> .env.ecr

.PHONY: local-ecr-login
local-ecr-login:
	if PROFILE=$$(aws configure list-profiles | grep single-cell-dev); then \
		aws ecr get-login-password --region us-west-2 --profile single-cell-dev | docker login --username AWS --password-stdin $$(aws sts get-caller-identity --profile single-cell-dev | jq -r .Account).dkr.ecr.us-west-2.amazonaws.com; \
	fi

.PHONY: local-init-test-data
local-init-test-data:
	docker-compose $(COMPOSE_OPTS) run --rm -T backend /bin/bash -c "pip3 install awscli && cd /single-cell-data-portal && scripts/setup_dev_data.sh"

.PHONY: local-init-host
local-init-host: oauth/pkcs12/certificate.pfx .env.ecr local-ecr-login local-start

.PHONY:
local-init: local-init-host local-init-test-data ## Launch a new local dev env and populate it with test data.

.PHONY: local-status
local-status: ## Show the status of the containers in the dev environment.
	docker ps -a | grep --color=no -e 'CONTAINER\|single-cell-data-portal'

.PHONY: local-rebuild
local-rebuild: .env.ecr local-ecr-login ## Rebuild local dev without re-importing data
	docker-compose $(COMPOSE_OPTS) build frontend backend processing wmg_processing database oidc localstack
	docker-compose $(COMPOSE_OPTS) up -d frontend backend processing database oidc localstack

local-rebuild-backend: .env.ecr local-ecr-login
	docker-compose $(COMPOSE_OPTS) build backend

local-rebuild-processing: .env.ecr local-ecr-login
	docker-compose $(COMPOSE_OPTS) build processing

local-rebuild-wmg-processing: .env.ecr local-ecr-login
	docker-compose $(COMPOSE_OPTS) build wmg_processing

local-rebuild-cellguide-pipeline: .env.ecr local-ecr-login
	docker-compose $(COMPOSE_OPTS) build cellguide_pipeline

.PHONY: local-sync
local-sync: local-rebuild local-init  ## Re-sync the local-environment state after modifying library deps or docker configs

.PHONY: local-start
local-start: .env.ecr ## Start a local dev environment that's been stopped.
	docker-compose $(COMPOSE_OPTS) up -d

.PHONY: local-stop
local-stop: ## Stop the local dev environment.
	docker-compose stop

.PHONY: local-clean
local-clean: ## Remove everything related to the local dev environment (including db data!)
	-if [ -f ./oauth/pkcs12/server.crt ] ; then \
	    export CERT=$$(docker run -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint "" soluto/oidc-server-mock:0.3.0 bash -c "openssl x509 -in server.crt -outform DER | sha1sum | cut -d ' ' -f 1"); \
	    echo ""; \
	    echo "Removing this certificate requires sudo access"; \
	    sudo security delete-certificate -Z $${CERT} /Library/Keychains/System.keychain; \
	fi;
	-rm -rf ./oauth/pkcs12/server*
	-rm -rf ./oauth/pkcs12/certificate*
	docker-compose rm -sf
	-docker volume rm single-cell-data-portal_database
	-docker volume rm single-cell-data-portal_localstack
	-docker network rm single-cell-data-portal_corporanet
	-docker network rm single-cell-data-portal_default

.PHONY: local-logs
local-logs: ## Tail the logs of the dev env containers. ex: make local-logs CONTAINER=backend
	docker-compose logs -f $(CONTAINER)

.PHONY: local-shell
local-shell: ## Open a command shell in one of the dev containers. ex: make local-shell CONTAINER=frontend
	docker-compose exec $(CONTAINER) bash

.PHONY: local-unit-test
local-unit-test: local-unit-test-backend local-unit-test-wmg-backend local-unit-test-wmg-processing local-unit-test-cellguide-pipeline local-unit-test-processing local-unit-test-cxg-admin
# Run all backend and processing unit tests in the dev environment, with code coverage

.PHONY: local-unit-test-backend
local-unit-test-backend: 
	docker-compose run --rm -T backend bash -c \
	"cd /single-cell-data-portal && coverage run  $(COVERAGE_RUN_ARGS) -m pytest --alluredir=./allure-results tests/unit/backend/layers/ tests/unit/backend/common/";

.PHONY: local-unit-test-wmg-backend
local-unit-test-wmg-backend: 
	docker-compose run --rm -T backend bash -c \
	"cd /single-cell-data-portal && coverage run $(COVERAGE_RUN_ARGS) -m pytest --alluredir=./allure-results tests/unit/backend/wmg/";

.PHONY: local-integration-test-backend
local-integration-test-backend:
	docker-compose run --rm -e INTEGRATION_TEST=true -e DB_URI=postgresql://corpora:test_pw@database -T backend \
	bash -c "cd /single-cell-data-portal && coverage run $(COVERAGE_RUN_ARGS) -m pytest tests/unit/backend/layers/ tests/unit/backend/common/";

.PHONY: local-unit-test-processing
local-unit-test-processing: # Run processing-unittest target in `processing` Docker container
	docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -T processing \
	bash -c "cd /single-cell-data-portal && coverage run $(COVERAGE_RUN_ARGS) -m pytest --alluredir=./allure-results tests/unit/processing/";

.PHONY: local-unit-test-wmg-processing
local-unit-test-wmg-processing: # Run processing-unittest target in `wmg_processing` Docker container
	echo "Running all wmg processing unit tests"; \
	docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -T wmg_processing \
	bash -c "cd /single-cell-data-portal && make wmg-processing-unittest;"

.PHONY: local-unit-test-cellguide-pipeline
local-unit-test-cellguide-pipeline: # Run processing-unittest target in `cellguide_pipeline` Docker container
	echo "Running all cellguide pipeline unit tests"; \
	docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -T cellguide_pipeline \
	bash -c "cd /single-cell-data-portal && make cellguide-pipeline-unittest;"	

.PHONY: local-unit-test-cxg-admin
local-unit-test-cxg-admin:
	docker-compose run --rm -T backend bash -c \
	"cd /single-cell-data-portal && coverage run  $(COVERAGE_RUN_ARGS) -m pytest --alluredir=./allure-results tests/unit/scripts/";

# We optionally pass BOTO_ENDPOINT_URL if it is set, even if it is
# set to be the empty string.
# Note that there is a distinction between BOTO_ENDPOINT_URL being
# the empty string (in which case we override the existing variable
# defined in docker-compose.yml to be empty string), and not being
# set (in which case the default from docker-compose is untouched)
.PHONY: local-functional-test
local-functional-test: ## Run functional tests in the dev environment
	if [ -n "$${BOTO_ENDPOINT_URL+set}" ]; then \
		EXTRA_ARGS="-e BOTO_ENDPOINT_URL"; \
	fi; \
	chamber -b secretsmanager exec corpora/backend/$${DEPLOYMENT_STAGE}/auth0-secret -- \
		docker-compose $(COMPOSE_OPTS) run --rm -T \
		-e CLIENT_ID -e CLIENT_SECRET \
		-e FUNCTEST_ACCOUNT_USERNAME -e FUNCTEST_ACCOUNT_PASSWORD \
		-e TEST_AUTH0_USER_ACCOUNT_PASSWORD -e TEST_APP_ID -e TEST_APP_SECRET \
		-e SUPER_CURATOR_API_KEY -e AUTH0_DOMAIN \
		-e DEPLOYMENT_STAGE  -e STACK_NAME -e API_BASE_URL $${EXTRA_ARGS} \
		backend bash -c "cd /single-cell-data-portal && make container-functionaltest"

.PHONY: local-smoke-test
local-smoke-test: ## Run frontend/e2e tests in the dev environment
	docker-compose $(COMPOSE_OPTS) run --rm -T frontend make smoke-test-with-local-dev


.PHONY: local-dbconsole
local-dbconsole: ## Connect to the local postgres database.
	psql "postgresql://corpora:test_pw@localhost:5432"

.PHONY: local-uploadjob
local-uploadjob: .env.ecr ## Run the upload task with a dataset_id and dropbox_url
	docker-compose $(COMPOSE_OPTS) run --rm -T -e DATASET_ID=$(DATASET_ID) -e DROPBOX_URL=$(DROPBOX_URL) processing sh -c "rm -rf /local.* && python3 -m backend.corpora.dataset_processing.process"

.PHONY: local-uploadfailure
local-uploadfailure: .env.ecr ## Run the upload failure lambda with a dataset id and cause
	docker-compose $(COMPOSE_OPTS) up -d upload_failures
	curl -v -XPOST "http://127.0.0.1:9000/2015-03-31/functions/function/invocations" -d '{"dataset_id": "$(DATASET_ID)", "error": {"Cause": "$(CAUSE)"}}'

.PHONY: local-uploadsuccess
local-uploadsuccess: .env.ecr ## Run the upload success lambda with a dataset id and cause
	docker-compose $(COMPOSE_OPTS) up -d upload_success
	curl -v -XPOST "http://127.0.0.1:9001/2015-03-31/functions/function/invocations" -d '{"dataset_id": "$(DATASET_ID)"}'

.PHONY: local-cxguser-cookie
local-cxguser-cookie: ## Get cxguser-cookie
	docker-compose $(COMPOSE_OPTS) run --rm backend bash -c "cd /single-cell-data-portal && python login.py"

.PHONY: coverage/combine
coverage/combine:
	- docker-compose $(COMPOSE_OPTS) run --rm -T backend bash -c "cd /single-cell-data-portal && coverage combine --data-file=$(COVERAGE_DATA_FILE)"

.PHONY: coverage/report
coverage/report-xml: coverage/combine
	docker-compose $(COMPOSE_OPTS) run --rm -T backend bash -c "cd /single-cell-data-portal && coverage xml --data-file=$(COVERAGE_DATA_FILE) -i --skip-empty"

.PHONY: coverage/report
coverage/report-html: coverage/combine
	docker-compose $(COMPOSE_OPTS) run --rm -T backend bash -c "cd /single-cell-data-portal && coverage html --data-file=$(COVERAGE_DATA_FILE) -i --skip-empty"
