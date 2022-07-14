variable image {
  type        = string
  description = "Image name"
}

variable task_role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
}

variable execution_role {
  type        = string
  description = "Execution role to use for fargate tasks - required for fargate services!"
  default     = ""
}

variable cmd {
  type        = list(string)
  description = "Command to run"
  default     = []
}

variable custom_stack_name {
  type        = string
  description = "Please provide the stack name"
}

variable remote_dev_prefix {
  type        = string
  description = "S3 storage path / db schema prefix"
  default     = ""
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
}
