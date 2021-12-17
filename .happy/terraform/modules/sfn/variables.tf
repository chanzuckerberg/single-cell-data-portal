variable job_definition_arn {
  type        = string
  description = "ARN of the batch job definition"
}

variable job_queue_arn {
  type        = string
  description = "ARN of the batch job queue"
}

variable lambda_success_handler {
  type        = string
  description = "ARN for the Lambda success handler"
}

variable lambda_error_handler {
  type        = string
  description = "ARN for the Lambda error handler"
}

variable role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
}

variable custom_stack_name {
  type        = string
  description = "Please provide the stack name"
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
  default     = "test"
}
