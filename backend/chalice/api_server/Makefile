SHELL=/bin/bash -o pipefail

# This is set for build reproducibility
export ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
export APP_NAME=corpora-api
export EXPORT_ENV_VARS_TO_LAMBDA=APP_NAME DEPLOYMENT_STAGE
export PYTHONHASHSEED=0
export S3_DEPLOYMENT_FILE=s3://org-dcp-infra-$(ACCOUNT_ID)/chalice/envs/$(DEPLOYMENT_STAGE)/deployed.json
export LOCAL_DEPLOYED_PATH=.chalice/deployed
export LOCAL_DEPLOYED_FILE=$(LOCAL_DEPLOYED_PATH)/$(DEPLOYMENT_STAGE).json
ifndef DEPLOYMENT_STAGE
$(error Please set DEPLOYMENT_STAGE in environment before running make commands)
endif

.PHONY: clean
clean:
	rm -rf dist .chalice/deployments chalicelib

.PHONY: destroy
destroy:
	chalice delete --stage $(DEPLOYMENT_STAGE)
	make remove-deployed


.PHONY: deploy
deploy: package
	chalice deploy --stage $(DEPLOYMENT_STAGE) --api-gateway-stage $(DEPLOYMENT_STAGE)
	make upload-deployed

.PHONY: ci-deploy
ci-deploy: retrieve-deployed
	make deploy

.PHONY: package
package: build-chalice-config
	mkdir -p chalicelib/code
	cp -R ../../corpora/common chalicelib/corpora
	cp -a $(APP_NAME).yml chalicelib
	rm -rf vendor dist/deployment
	mkdir vendor
	find vendor -name '*.pyc' -delete
	shopt -s nullglob; for wheel in vendor.in/*/*.whl; do pip install --target vendor --upgrade $$wheel; done

.PHONY: build-chalice-config
build-chalice-config:
	envsubst < ../../config/iam-policy-templates/$(APP_NAME)-lambda.json > .chalice/policy-$(DEPLOYMENT_STAGE).json
	cd .chalice; envsubst < config.in.json > config.json
	cd .chalice; for var in $$EXPORT_ENV_VARS_TO_LAMBDA; do \
            jq .stages.$(DEPLOYMENT_STAGE).environment_variables.$$var=env.$$var config.json | sponge config.json; \
        done
	cd .chalice; V=$$(git describe --tags --always) jq .stages.$(DEPLOYMENT_STAGE).environment_variables.VERSION=env.V config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).iam_policy_file=\"policy-\"+env.DEPLOYMENT_STAGE+\".json\" config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).tags.env=env.DEPLOYMENT_STAGE config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).tags.service=env.APP_NAME config.json | sponge config.json
	cd .chalice; jq .stages.$(DEPLOYMENT_STAGE).api_gateway_stage=env.DEPLOYMENT_STAGE config.json | sponge config.json

.PHONY: local-server
local-server: package
	python local_server.py

.PHONY: retrieve-deployed
retrieve-deployed:
	if [[ ! -d "$(LOCAL_DEPLOYED_PATH)" ]]; then mkdir $(LOCAL_DEPLOYED_PATH); fi
	aws s3 cp $(S3_DEPLOYMENT_FILE) $(LOCAL_DEPLOYED_FILE)

.PHONY: upload-deployed
upload-deployed:
	aws s3 cp $(LOCAL_DEPLOYED_FILE) $(S3_DEPLOYMENT_FILE)

.PHONY: remove-deployed
remove-deployed:
	aws s3 rm $(S3_DEPLOYMENT_FILE)
