terraform {
  required_version = "=0.12.20"

  backend "s3" {
    key     = "terraform/envs/dev/browser_backend.tfstate"
    encrypt = true
    region  = "us-east-1"
    profile = "czi-hca-dev"
  }
}

provider "aws" {
  version = "~> 2.33"
  region  = "us-east-1"
  profile = "czi-hca-dev"
}

provider "external" {
  version = "~> 1.1"
}

provider "template" {
  version = "~> 2.1"
}

module "browser" {
  source = "./terraform"
}
