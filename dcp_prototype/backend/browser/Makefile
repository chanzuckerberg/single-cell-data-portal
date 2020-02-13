SHELL=/bin/bash -o pipefail
tests:=$(wildcard tests/test*.py)

ifndef STAGE
$(error Please run "source environment" in the repo root directory before running make commands)
endif

deploy: init package
	export LAMBDA_SHA=$(shell sha256sum dist/deployment.zip | cut -f 1 -d ' ' | base64)
	aws s3 cp dist/deployment.zip s3://$(TF_S3_BUCKET)/$(LAMBDA_SHA).zip
	jq -f scripts/preprocess_chalice_tf_config.jq dist/chalice.tf.json > terraform/chalice.tf.json
	terraform apply

plan: init package
	$(eval export LAMBDA_SHA=$(shell sha256sum dist/deployment.zip | cut -f 1 -d ' ' | base64))
	aws s3 cp dist/deployment.zip s3://$(TF_S3_BUCKET)/$(LAMBDA_SHA).zip
	jq -f scripts/preprocess_chalice_tf_config.jq dist/chalice.tf.json > terraform/chalice.tf.json
	terraform plan

# init-tf prepares the repo for Terraform commands. It assembles the partial S3 backend config as a JSON file, `aws_config.json`.
# This file is referenced by the TF_CLI_ARGS_init environment variable, which is set by running `source environment`.
init:
	-rm -f $(TF_DATA_DIR)/*.tfstate
	mkdir -p $(TF_DATA_DIR)
	jq -n ".region=env.AWS_REGION | .bucket=env.TF_S3_BUCKET | .key=env.APP_NAME+env.STAGE" > $(TF_DATA_DIR)/aws_config.json
	terraform init

build-chalice-config:
	envsubst < iam/policy-templates/$(APP_NAME)-lambda.json > .chalice/policy-$(STAGE).json
	cd .chalice; jq .app_name=env.APP_NAME < config.in.json > config.json
	cd .chalice; for var in $$EXPORT_ENV_VARS_TO_LAMBDA; do \
            jq .stages.$(STAGE).environment_variables.$$var=env.$$var config.json | sponge config.json; \
        done
	cd .chalice; V=$$(git describe --tags --always) jq .stages.$(STAGE).environment_variables.VERSION=env.V config.json | sponge config.json
	cd .chalice; jq .stages.$(STAGE).tags.env=env.STAGE config.json | sponge config.json
	cd .chalice; jq .stages.$(STAGE).tags.service=env.APP_NAME config.json | sponge config.json
	cd .chalice; jq .stages.$(STAGE).tags.owner=env.OWNER config.json | sponge config.json
	cd .chalice; jq .stages.$(STAGE).api_gateway_stage=env.STAGE config.json | sponge config.json

# package prepares a Lambda zipfile with the help of the Chalice packager (which also emits a SAM template).
# We also inject any wheels found in vendor.in, and rewrite the zipfile to make the build reproducible.
package: build-chalice-config
	rm -rf vendor dist/deployment
	mkdir vendor
	cp -a $(APP_NAME) $(APP_NAME)-api.yml vendor
	find vendor -name '*.pyc' -delete
	find vendor -exec touch -t 201901010000 {} \; # Reset mtimes on all vendor files to make zipfile contents reproducible
	shopt -s nullglob; for wheel in vendor.in/*/*.whl; do pip install --target vendor --upgrade $$wheel; done
	chalice package --pkg-format terraform --stage $(STAGE) dist
	cd dist; mkdir deployment; cd deployment; unzip -q -o ../deployment.zip
	$(MAKE) prune
	find dist/deployment -exec touch -t 201901010000 {} \; # Reset mtimes on all dep files to make zipfile contents reproducible
	rm dist/deployment.zip
	cd dist/deployment; zip -q -X -r ../deployment.zip .

# Remove unused dependencies that consume substantial space in the package.
# See https://docs.aws.amazon.com/lambda/latest/dg/lambda-runtimes.html
prune:
	cd dist/deployment; rm -rf awscli* boto3* botocore* cryptography* swagger_ui_bundle/vendor/swagger-ui-2* connexion/vendor/swagger-ui*

destroy: init
	terraform destroy

clean:
	rm -rf dist .terraform* .chalice/deployments

test: package lint $(tests)
	coverage combine
	rm -f .coverage.*

$(tests): %.py : lint
	coverage run -p --source=browser $*.py -v

unittest: package
	python -m unittest discover -p test_*api.py -t . -v

.PHONY: deploy package init destroy clean lint test
