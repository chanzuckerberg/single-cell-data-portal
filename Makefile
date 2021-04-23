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


.PHONY: fmt
fmt:
	black .

.PHONY: lint
lint:
	flake8 backend tests

.PHONY: unit-test
unit-test: local-unit-test
	# Keeping old target name for reverse comatibility

.PHONY: container-unittest
container-unittest:
	# This target is intended to be run INSIDE a container
	DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning python3 -m coverage run \
		-m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose;

.PHONY: processing-unittest
processing-unittest:
	# This target is intended to be run INSIDE a container
	DEPLOYMENT_STAGE=test PYTHONWARNINGS=ignore:ResourceWarning python3 -m coverage run \
		-m unittest discover --start-directory tests/unit/processing_container --top-level-directory . --verbose;

.PHONY: functional-test
functional-test: local-functional-test
	# Keeping old target name for reverse comatibility

.PHONY: container-functionaltest
container-functionaltest:
	# This target is intended to be run INSIDE a container
	python3 -m unittest discover --start-directory tests/functional --top-level-directory . --verbose

.PHONY: local-backend
local-backend:
	$(MAKE) local-server -C ./backend/chalice/api_server

.PHONY: smoke-test-prod-build
smoke-test-prod-build:
	$(MAKE) smoke-test-prod-build -C ./frontend

.PHONY: smoke-test-with-local-backend
smoke-test-with-local-backend:
	$(MAKE) smoke-test-with-local-backend -C ./frontend

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

.env.ecr:
	echo DOCKER_REPO=$$(aws sts get-caller-identity --profile single-cell-dev | jq -r .Account).dkr.ecr.us-west-2.amazonaws.com/ > .env.ecr;

.PHONY: local-ecr-login
local-ecr-login:
	if PROFILE=$$(aws configure list-profiles | grep single-cell-dev); then \
		aws ecr get-login-password --region us-west-2 --profile single-cell-dev | docker login --username AWS --password-stdin $$(aws sts get-caller-identity --profile single-cell-dev | jq -r .Account).dkr.ecr.us-west-2.amazonaws.com; \
	fi

.PHONY: local-hostconfig
local-hostconfig:
	if [ "$$(uname -s)" == "Darwin" ]; then \
	  sudo ./scripts/happy hosts install; \
	fi

.PHONY: local-nohostconfig
local-nohostconfig:
	if [ "$$(uname -s)" == "Darwin" ]; then \
	  sudo ./scripts/happy hosts uninstall; \
	fi

.PHONY: local-init
local-init: oauth/pkcs12/certificate.pfx .env.ecr local-ecr-login local-hostconfig ## Launch a new local dev env and populate it with test data.
	docker-compose $(COMPOSE_OPTS) up -d frontend backend database oidc localstack
	docker-compose $(COMPOSE_OPTS) run --rm -T backend /bin/bash -c "pip3 install awscli && /corpora-data-portal/scripts/setup_dev_data.sh"

.PHONY: local-status
local-status: ## Show the status of the containers in the dev environment.
	docker ps -a | grep --color=no -e 'CONTAINER\|corpora-data-portal'

.PHONY: local-rebuild
local-rebuild: .env.ecr local-ecr-login ## Rebuild local dev without re-importing data
	docker-compose $(COMPOSE_OPTS) build frontend backend
	docker-compose $(COMPOSE_OPTS) up -d frontend backend database oidc localstack

.PHONY: local-sync
local-sync: local-rebuild local-init local-hostconfig ## Re-sync the local-environment state after modifying library deps or docker configs

.PHONY: local-start
local-start: .env.ecr ## Start a local dev environment that's been stopped.
	docker-compose $(COMPOSE_OPTS) up -d

.PHONY: local-stop
local-stop: ## Stop the local dev environment.
	docker-compose stop

.PHONY: local-clean
local-clean: local-nohostconfig ## Remove everything related to the local dev environment (including db data!)
	-if [ -f ./oauth/pkcs12/server.crt ] ; then \
	    export CERT=$$(docker run -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint "" soluto/oidc-server-mock:0.3.0 bash -c "openssl x509 -in server.crt -outform DER | sha1sum | cut -d ' ' -f 1"); \
	    echo ""; \
	    echo "Removing this certificate requires sudo access"; \
	    sudo security delete-certificate -Z $${CERT} /Library/Keychains/System.keychain; \
	fi;
	-rm -rf ./oauth/pkcs12/server*
	-rm -rf ./oauth/pkcs12/certificate*
	docker-compose rm -sf
	-docker volume rm corpora-data-portal_database
	-docker volume rm corpora-data-portal_localstack
	-docker network rm corpora-data-portal_corporanet
	-docker network rm corpora-data-portal_default

.PHONY: local-logs
local-logs: ## Tail the logs of the dev env containers. ex: make local-logs CONTAINER=backend
	docker-compose logs -f $(CONTAINER)

.PHONY: local-shell
local-shell: ## Open a command shell in one of the dev containers. ex: make local-shell CONTAINER=frontend
	docker-compose exec $(CONTAINER) bash

.PHONY: local-unit-test
local-unit-test: ## Run backend tests in the dev environment
	@if [ -z "$(path)" ]; then \
        echo "Running all tests"; \
		export CI=""; \
		ci_env=""; \
		if [ ! -z "$(CODECOV_TOKEN)" ]; then \
			ci_env=$$(bash <(curl -s https://codecov.io/env)); \
			CI=true; \
		fi; \
		$(MAKE) _unit_test; \
	else \
		echo "Running test(s): $(path)"; \
		docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -T backend bash -c "cd /corpora-data-portal && python -m unittest $(path)"; \
	fi

.PHONY: local-processing-unit-test
local-processing-unit-test: ## Run backend tests in the dev environment
	@if [ -z "$(path)" ]; then \
        echo "Running all tests"; \
		export CI=""; \
		ci_env=""; \
		if [ ! -z "$(CODECOV_TOKEN)" ]; then \
			ci_env=$$(bash <(curl -s https://codecov.io/env)); \
			CI=true; \
		fi; \
		$(MAKE) _processing_unit_test; \
	else \
		echo "Running test(s): $(path)"; \
		docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -T backend bash -c "cd /corpora-data-portal && python -m unittest $(path)"; \
	fi

.PHONY: _unit_test
_unit_test:
	docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -e CI $$ci_env -T backend bash -c "cd /corpora-data-portal && make container-unittest && if [ "$${CI}" = "true" ]; then apt-get update && apt-get install -y git && bash <(curl -s https://codecov.io/bash) -cF backend,python,unitTest; fi";

.PHONY: _processing_unit_test
_processing_unit_test:
	docker-compose $(COMPOSE_OPTS) run --rm -e DEV_MODE_COOKIES= -e CI $$ci_env -T processing bash -c "cd /corpora-data-portal && make processing-unittest && if [ "$${CI}" = "true" ]; then apt-get update && apt-get install -y git && bash <(curl -s https://codecov.io/bash) -cF backend,python,unitTest; fi";

# We optionally pass BOTO_ENDPOINT_URL if it is set, even if it is
# set to be the empty string.
# Note that there is a distinction between BOTO_ENDPOINT_URL being
# the empty string (in which case we override the existing variable
# defined in docker-compose.yml to be empty string), and not being
# set (in which case the default from docker-compose is untouched)
#
# Unfortunately, this isn't working properly if DEPLOYMENT_STAGE is not test.
# If you want to run this locally against staging, use the following mouthful of a
# command, replacing IAM_ROLE with a role you have access to:
# export IAM_ROLE=insert_iam_arn_here
# eval $(echo $(AWS_PROFILE=single-cell-dev aws sts assume-role --role-arn ${IAM_ROLE} --role-session-name functional-test \
#     | jq -r '"AWS_ACCESS_KEY_ID=\(.Credentials.AccessKeyId) AWS_SECRET_ACCESS_KEY=\(.Credentials.SecretAccessKey) AWS_SESSION_TOKEN=\(.Credentials.SessionToken)"') \
# 	AWS_REGION=us-west-2 BOTO_ENDPOINT_URL= DEPLOYMENT_STAGE=staging make local-functional-test)
.PHONY: local-functional-test
local-functional-test: export AWS_PROFILE=$(TEST_AWS_PROFILE)
local-functional-test: ## Run functional tests in the dev environment
	if [ -n "$${BOTO_ENDPOINT_URL+set}" ]; then \
		EXTRA_ARGS="-e BOTO_ENDPOINT_URL"; \
	fi; \
	chamber -b secretsmanager exec corpora/backend/$${DEPLOYMENT_STAGE}/auth0-secret -- \
		docker-compose $(COMPOSE_OPTS) run --rm -T -e CLIENT_ID -e CLIENT_SECRET -e TEST_ACCOUNT_USERNAME -e TEST_ACCOUNT_PASSWORD -e DEPLOYMENT_STAGE -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_SESSION_TOKEN $${EXTRA_ARGS} \
		backend bash -c "cd /corpora-data-portal && make container-functionaltest"

.PHONY: local-smoke-test
local-smoke-test: ## Run frontend/e2e tests in the dev environment
	docker-compose $(COMPOSE_OPTS) run --rm -T frontend make smoke-test-with-local-dev

.PHONY: local-e2e
local-e2e: ## Run frontend/e2e tests
	if [ -n "$${BOTO_ENDPOINT_URL+set}" ]; then \
		EXTRA_ARGS="-e BOTO_ENDPOINT_URL"; \
	fi; \
	docker-compose $(COMPOSE_OPTS) run --no-deps -e DEPLOYMENT_STAGE -e AWS_REGION -e AWS_DEFAULT_REGION -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_SESSION_TOKEN $${EXTRA_ARGS} -T frontend make e2e

.PHONY: local-dbconsole
local-dbconsole: ## Connect to the local postgres database.
	psql "postgresql://corpora:test_pw@localhost:5432"

.PHONY: local-uploadjob
local-uploadjob: .env.ecr ## Run the upload task with a dataset_id and dropbox_url
	docker-compose $(COMPOSE_OPTS) run --rm -T -e DATASET_ID=$(DATASET_ID) -e DROPBOX_URL=$(DROPBOX_URL) processing sh -c "rm -rf /local.* && python3 -m backend.corpora.dataset_processing.process"

.PHONY: local-uploadfailure
local-uploadfailure: .env.ecr ## Run the upload failure lambda with a dataset id and cause
	docker-compose $(COMPOSE_OPTS) up -d upload_failures
	curl -v -XPOST "http://127.0.0.1:9000/2015-03-31/functions/function/invocations" -d '{"dataset_uuid": "$(DATASET_UUID)", "error": {"Cause": "$(CAUSE)"}}'

.PHONY: local-cxguser-cookie
local-cxguser-cookie: ## Get cxguser-cookie
	docker-compose $(COMPOSE_OPTS) run --rm backend bash -c "cd /corpora-data-portal && python login.py"
