resource "aws_api_gateway_base_path_mapping" "mapping" {
  api_id      = var.api_gateway_id
  stage_name  = var.deployment_stage
  domain_name = var.cert_domain_name
}

resource "aws_api_gateway_domain_name" "api_domain" {
  domain_name              = var.cert_domain_name
  regional_certificate_arn = var.aws_acm_cert_arn

  endpoint_configuration {
    types = ["REGIONAL"]
  }
}

resource "aws_route53_record" "api_domain" {
  zone_id = var.aws_route53_zone_id
  name    = var.cert_domain_name
  type    = "CNAME"
  ttl     = "300"
  records = ["${aws_api_gateway_domain_name.api_domain.regional_domain_name}"]
}
