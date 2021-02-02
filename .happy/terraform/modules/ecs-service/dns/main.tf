# This template creates a route53 cname for a shared alb resource.
# 

data aws_route53_zone dns_record {
  name = var.zone
}

resource aws_route53_record dns_record_0 {
  name    = "${var.custom_stack_name}-${var.app_name}.${var.zone}"
  type    = "A"
  zone_id = data.aws_route53_zone.dns_record.zone_id
  alias {
    name                   = var.alb_dns
    zone_id                = var.canonical_hosted_zone
    evaluate_target_health = false
  }
}
