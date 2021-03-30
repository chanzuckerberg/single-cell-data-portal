SHELL:=/bin/bash
ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
BROWSER_S3=corpora-static-site-$(DEPLOYMENT_STAGE)-$(ACCOUNT_ID)
S3_ENVIRONMENT_FILE=s3://$(BROWSER_S3)/env.production

check-var-defined = $(if $(strip $($1)),,$(error "$1" is not defined))

deploy: init retrieve-vars

	npm run build
	aws s3 sync ./public s3://$(BROWSER_S3)/

.PHONY: init
init:
	npm install

.PHONY: lint
lint:
	npm run prettier-check
	npm run lint

.PHONY: clean
clean:
	npm run clean

.PHONY: retrieve-vars
retrieve-vars:
	$(call check-var-defined,DEPLOYMENT_STAGE)
	aws s3 cp $(S3_ENVIRONMENT_FILE) .env.production
	# Copy the appropriate configs
	cp src/configs/$(DEPLOYMENT_STAGE).js src/configs/configs.js

.PHONY: upload-vars
upload-vars:
	aws s3 cp .env.production $(S3_ENVIRONMENT_FILE)

.PHONY: smoke-test-prod-build
smoke-test-prod-build:
	npm run smoke-test-prod-build

.PHONY: smoke-test-with-local-backend
smoke-test-with-local-backend:
	npm run smoke-test-with-local-backend

.PHONY: e2e
e2e:
	TEST_ENV=$(DEPLOYMENT_STAGE) npm run e2e

.PHONY: smoke-test-with-local-dev
smoke-test-with-local-dev:
	npm run e2e; \
	test_result=$$?; \
	cp src/configs/local.js src/configs/configs.js; \
	exit $$test_result

