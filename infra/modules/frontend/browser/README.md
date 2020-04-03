**copied from [https://github.com/chanzuckerberg/cztack/tree/v0.30.0/aws-single-page-static-site](https://github.com/chanzuckerberg/cztack/tree/v0.30.0/aws-single-page-static-site)**

# AWS S3 Single Page Static Site

This module will create a public static site powered by [S3](https://aws.amazon.com/s3/) with [Cloudfront](https://aws.amazon.com/cloudfront/) in front.

A few things this module does for you–

* creates a private s3 bucket
  * ensures data is encrypted by default using S3 Server Side Encryption
  * restricts access to a cloudfront [origin access identity](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/private-content-restricting-access-to-s3.html)
  identity
* creates a cloudfront distribution
  * TLS required through [ACM](https://aws.amazon.com/certificate-manager/)
  * allows you to configure multiple aliases
  * configures the distribution to serve a single page application
  * configures gzip compression
* configures route53 records for your distribution
* tags all resources with our [consistent tags](../README.md#Consistent%20Tagging)

This module is specifically designed for Single Page Applications– In order to make this work, cloudfront redirects certain S3 http responses (`403`, `404`, `503`) back to your index page for handling.

## Example

```hcl
module "site" {
  source = "github.com/chanzuckerberg/cztack//aws-single-page-static-site?ref=v0.14.0"

  cert_domain_name               = "..."
  cert_subject_alternative_names = "..."
  aws_route53_zone_id            = "..."
  validation_record_ttl          = "..."

  # Variables used for tagging
  env     = "..."
  project = "..."
  service = "..."
  owner   = "..."
}
```

<!-- START -->
## Providers

| Name | Version |
|------|---------|
| aws | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:-----:|
| aliases | Vanity aliases. Make sure your provided cert supports these. | `list` | `[]` | no |
| aws\_acm\_cert\_arn | An AWS ACM cert. Note that Cloudfront requires certs to be in us-east-1. | `string` | n/a | yes |
| aws\_route53\_zone\_id | A route53 zone ID used to write records. | `string` | n/a | yes |
| bucket\_name | Name of the bucket to created. If not given, it will use the domain name. | `string` | `""` | no |
| cloudfront\_price\_class | Cloudfront [price class](https://aws.amazon.com/cloudfront/pricing/). | `string` | `"PriceClass_100"` | no |
| env | Env for tagging and naming. See [doc](../README.md#consistent-tagging) | `string` | n/a | yes |
| index\_document\_path | The path to the index document of your site. | `string` | `"index.html"` | no |
| minimum\_tls\_version | Minimum TLS version to accept. | `string` | `"TLSv1.1_2016"` | no |
| owner | Owner for tagging and naming. See [doc](../README.md#consistent-tagging) | `string` | n/a | yes |
| path\_pattern | The pattern (for example, images/\*.jpg) that specifies which requests you want this cache behavior to apply to. | `string` | `"*"` | no |
| project | Project for tagging and naming. See [doc](../README.md#consistent-tagging) | `string` | n/a | yes |
| service | Service for tagging and naming. See [doc](../README.md#consistent-tagging) | `string` | n/a | yes |
| subdomain | The subdomain for this static site. | `string` | n/a | yes |
| refer_secret | A secret string passed from cloudfront to s3 to prevent public access to the bucket| `string` | n/a | yes|

## Outputs

| Name | Description |
|------|-------------|
| bucket\_arn | n/a |
| bucket\_name | n/a |
| cloudfront\_arn | n/a |
| cloudfront\_domain\_name | n/a |
| cloudfront\_hosted\_zone\_id | n/a |

<!-- END -->
