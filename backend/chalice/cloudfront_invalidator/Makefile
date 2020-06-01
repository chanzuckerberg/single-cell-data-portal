SHELL=/bin/bash -o pipefail

# This is set for build reproducibility
export ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
export APP_NAME?=browser-frontend
ifndef DEPLOYMENT_STAGE
$(error Please set DEPLOYMENT_STAGE in environment before running make commands)
endif
export EXPORT_ENV_VARS_TO_LAMBDA=APP_NAME DEPLOYMENT_STAGE

.PHONY: clean
clean:
	rm -rf dist .chalice/deployments chalicelib .chalice/config.json .chalice/policy-*

.PHONY: destroy
destroy:
	chalice delete --stage $(DEPLOYMENT_STAGE)

.PHONY: package
package: build-chalice-config
	mkdir -p chalicelib/code
	cp -R ../../code chalicelib

.PHONY: deploy
deploy: package
	chalice deploy --stage $(DEPLOYMENT_STAGE)

.PHONY: build-chalice-config
build-chalice-config:
	envsubst < ../../iam-policy-templates/$(APP_NAME)-lambda.json > .chalice/policy-$(DEPLOYMENT_STAGE).json
	cd .chalice; envsubst < config.in.json > config.json
	cd .chalice; for var in $$EXPORT_ENV_VARS_TO_LAMBDA; do \
            jq .stages.$(DEPLOYMENT_STAGE).environment_variables.$$var=env.$$var config.json | sponge config.json; \
        done
	cd .chalice; V=$$(git describe --tags --always) jq .stages.$(DEPLOYMENT_STAGE).environment_variables.VERSION=env.V config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).iam_policy_file=\"policy-\"+env.DEPLOYMENT_STAGE+\".json\" config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).tags.env=env.DEPLOYMENT_STAGE config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).tags.service=env.APP_NAME config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).api_gateway_stage=env.DEPLOYMENT_STAGE config.json | sponge config.json
