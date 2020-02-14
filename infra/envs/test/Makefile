ACCOUNT_ID=$(shell aws sts get-caller-identity --query Account --output text)
DEPLOYMENT_STAGE=test
TF_STATE_BUCKET=org-dcp-infra-856863240400
S3_TFVARS_FILE=s3://$(TF_STATE_BUCKET)/terraform/envs/$(DEPLOYMENT_STAGE)/terraform.tfvars

default: plan

init:
	envsubst < main.tf.in > main.tf
	terraform init \
    -backend-config="bucket=org-dcp-infra-${ACCOUNT_ID}"

plan: retrieve-vars
	terraform plan -detailed-exitcode

apply: retrieve-vars
	terraform apply --backup=-

destroy: state/rm
	terraform destroy --backup=-

retrieve-vars:
	aws s3 cp $(S3_TFVARS_FILE) terraform.tfvars

upload-vars:
	aws s3 cp terraform.tfvars $(S3_TFVARS_FILE)
