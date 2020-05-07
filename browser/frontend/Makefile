SHELL:=/bin/bash
ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
BROWSER_S3=dcp-static-site-$(DEPLOYMENT_STAGE)-$(ACCOUNT_ID)
S3_ENVIRONMENT_FILE=s3://$(BROWSER_S3)/env.production

ifndef DEPLOYMENT_STAGE
$(error Please set the DEPLOYMENT_STAGE environment before deploying)
endif

deploy: init retrieve-vars
	npm run build
	aws s3 sync ./public s3://$(BROWSER_S3)/

.PHONY: init
init:
	npm install

.PHONY: clean
clean:
	npm run clean

.PHONY: retrieve-vars
retrieve-vars:
	aws s3 cp $(S3_ENVIRONMENT_FILE) .env.production

.PHONY: upload-vars
upload-vars:
	aws s3 cp .env.production $(S3_ENVIRONMENT_FILE)


.PHONY: development
develop: init
	aws s3 cp s3://$(BROWSER_S3)/local/env.development.local .env.development
	gatsby develop

.PHONY: development
upload-var-develop:
	aws s3 cp .env.development s3://$(BROWSER_S3)/local/env.development.local
