data "aws_caller_identity" "current" {}

output "account_id" {
  value = "${data.aws_caller_identity.current.account_id}"
}

terraform {
  required_version = "=0.12.20"

  backend "s3" {
    key     = "terraform/envs/${var.deployment_stage}/state.tfstate"
    encrypt = true
    region  = "us-east-1"
    profile = "czi-hca-dev"
  }
}

provider "aws" {
  version = "~> 2.0"
  region  = "us-east-1"
  profile = "czi-hca-dev"
}

module "ledger" {
  source = "../../modules/backend/ledger"

  deployment_stage = "${var.deployment_stage}"

  // Database
  db_username = "${var.db_username}"
  db_password = "${var.db_password}"
  db_instance_count = "${var.db_instance_count}"
  preferred_maintenance_window = "${var.preferred_maintenance_window}"
}