variable "project" {
  type        = string
  description = "Project for tagging and naming. See [doc](../README.md#consistent-tagging)"
}

variable "env" {
  type        = string
  description = "Env for tagging and naming. See [doc](../README.md#consistent-tagging)"
}

variable "service" {
  type        = string
  description = "Service for tagging and naming. See [doc](../README.md#consistent-tagging)"
}

variable "owner" {
  type        = string
  description = "Owner for tagging and naming. See [doc](../README.md#consistent-tagging)"
}

variable "subdomain" {
  type        = string
  description = "The subdomain for this static site."
}

variable "aws_acm_cert_arn" {
  type        = string
  description = "An AWS ACM cert. Note that Cloudfront requires certs to be in us-east-1. "
}

variable "index_document_path" {
  type        = string
  default     = "index.html"
  description = "The path to the index document of your site."
}

variable "aws_route53_zone_id" {
  type        = string
  description = "A route53 zone ID used to write records."
}

variable "cloudfront_price_class" {
  type        = string
  default     = "PriceClass_100"
  description = "Cloudfront [price class](https://aws.amazon.com/cloudfront/pricing/)."
}

variable "minimum_tls_version" {
  type        = string
  default     = "TLSv1.1_2016"
  description = "Minimum TLS version to accept."
}

variable "aliases" {
  type        = list
  default     = []
  description = "Vanity aliases. Make sure your provided cert supports these."
}

variable "bucket_name" {
  type        = string
  description = "Name of the bucket to created. If not given, it will use the domain name."
  default     = ""
}

variable "path_pattern" {
  type        = string
  description = "The pattern (for example, images/*.jpg) that specifies which requests you want this cache behavior to apply to."
  default     = "*"
}
variable "refer_secret" {
  type        = string
  description = "Used by cloudfront to restrict direct access to the s3 bucket serving the static website. Should be a sufficiently long string."
}
