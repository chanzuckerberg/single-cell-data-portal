output dns_prefix {
  value       = "${var.custom_stack_name}-${var.app_name}"
  description = "User-facing URL for this service."
}
