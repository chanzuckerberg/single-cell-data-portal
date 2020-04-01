data "aws_caller_identity" "current" {}

output "account_id" {
  value = "${data.aws_caller_identity.current.account_id}"
}

terraform {
  required_version = "=0.12.20"

  backend "s3" {
    encrypt = true
    region  = "us-east-1"
    profile = "single-cell-dev"
  }
}

provider "aws" {
  version = "~> 2.0"
  region  = "us-east-1"
  profile = "single-cell-dev"
}

//module "ledger" {
//  source = "../../modules/backend/ledger"
//
//  deployment_stage = "${var.deployment_stage}"
//
//  // Database
//  db_username                  = "${var.db_username}"
//  db_password                  = "${var.db_password}"
//  db_instance_count            = "${var.db_instance_count}"
//  preferred_maintenance_window = "${var.preferred_maintenance_window}"
//}

//module "browser" {
//  source = "../../modules/backend/browser"
//
//  deployment_stage = "${var.deployment_stage}"
//
//  // Database
//  db_username                  = "${var.browser_db_username}"
//  db_password                  = "${var.browser_db_password}"
//  db_instance_count            = "${var.browser_db_instance_count}"
//  preferred_maintenance_window = "${var.browser_preferred_maintenance_window}"
//}

module "browser_site_cert" {
  source = "github.com/chanzuckerberg/cztack//aws-acm-cert?ref=v0.29.0"

  # the cert domain name
  cert_domain_name               = "browser-testing.dev.single-cell.czi.technology"
  cert_subject_alternative_names = { "www.browser-testing.dev.single-cell.czi.technology" = "Z3GU3D81Z7R7CN" }

  # the route53 zone for validating the `cert_domain_name`
  aws_route53_zone_id = "Z3GU3D81Z7R7CN"

  # variables for tags
  env     = var.deployment_stage
  project = "single-cell"
  service = "browser"
  owner   = "czi-single-cell"
}

module "browser_frontend" {
  source = "../../modules/frontend/browser"

  aws_route53_zone_id = "Z3GU3D81Z7R7CN"
  aws_acm_cert_arn    = module.browser_site_cert.arn
  bucket_name         = "dcp-static-site-${var.deployment_stage}-${data.aws_caller_identity.current.account_id}"
  subdomain           = "browser-testing"

  # Variables used for tagging
  env     = var.deployment_stage
  project = "single-cell"
  service = "browser"
  owner   = "czi-single-cell"
  refer_secret = var.refer_secret
}

