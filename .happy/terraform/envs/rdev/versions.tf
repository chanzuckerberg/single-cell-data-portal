terraform {
  required_version = "~> 1.3.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 3.75.2"
    }
  }

  backend "s3" {
    bucket         = "terragrunt-engine-state"
    dynamodb_table = "terragrunt-engine-state-lock"
    key            = "terraform/single-cell-data-portal/envs/rdev/components/rdevstack.tfstate"
    encrypt        = true
    region         = "us-west-2"
    role_arn       = "arn:aws:iam::533267185808:role/tfe-si"
  }
}
