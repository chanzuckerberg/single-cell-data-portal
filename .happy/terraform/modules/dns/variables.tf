variable custom_stack_name {
  type        = string
  description = "Please provide the stack name"
}

variable app_name {
  type        = string
  description = "Please provide the ECS service name"
}

variable zone {
  type        = string
  description = "Route53 zone name. Trailing . must be OMITTED!"
}

variable alb_dns {
  type        = string
  description = "DNS name for the shared ALB"
}

variable canonical_hosted_zone {
  type        = string
  description = "Route53 zone for the shared ALB"
}
