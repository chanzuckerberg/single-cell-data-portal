ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
DEPLOYMENT_STAGE=staging
TF_STATE_BUCKET=org-dcp-infra-$(ACCOUNT_ID)
S3_TFVARS_PREFIX=s3://$(TF_STATE_BUCKET)/terraform/envs/$(DEPLOYMENT_STAGE)
S3_TFVARS_FILE=$(S3_TFVARS_PREFIX)/terraform.tfvars

default: plan

.PHONY: init
init:
	terraform init \
    -backend-config="bucket=org-dcp-infra-${ACCOUNT_ID}" \
    -backend-config="key=terraform/envs/$(DEPLOYMENT_STAGE)/state.tfstate"

.PHONY: plan
plan: retrieve-vars
	terraform plan -detailed-exitcode

.PHONY: apply
apply: retrieve-vars
	terraform apply --backup=-

.PHONY: destroy
destroy:
	terraform destroy --backup=-

.PHONY: rm-state
rm-state:
	aws s3 rm --recursive $(S3_TFVARS_PREFIX)

.PHONY: retrieve-vars
retrieve-vars:
	aws s3 cp $(S3_TFVARS_FILE) terraform.tfvars

.PHONY: upload-vars
upload-vars:
	aws s3 cp terraform.tfvars $(S3_TFVARS_FILE)
