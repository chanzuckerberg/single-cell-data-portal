ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
DEPLOYMENT_STAGE=dev
TF_STATE_BUCKET=org-dcp-infra-$(ACCOUNT_ID)
S3_TFVARS_FILE=s3://$(TF_STATE_BUCKET)/terraform/envs/$(DEPLOYMENT_STAGE)/terraform.tfvars

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
destroy: state/rm
	terraform destroy --backup=-

.PHONY: retrieve-vars
retrieve-vars:
	aws s3 cp $(S3_TFVARS_FILE) terraform.tfvars

.PHONY: upload-vars
upload-vars:
	aws s3 cp terraform.tfvars $(S3_TFVARS_FILE)
