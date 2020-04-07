variable "deployment_stage" {
  type = string
}

// RDS

variable "db_username" {
  type = string
}
variable "db_password" {
  type = string
}
variable "db_instance_count" {
  type = string
}
variable "preferred_maintenance_window" {
  type = string
}

// API Gateway
variable "aws_acm_cert_arn" {
  type= string
}
variable "cert_domain_name" {
  type = string
}
variable "aws_route53_zone_id" {
  type = string
}

