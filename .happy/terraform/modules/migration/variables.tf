variable image {
  type        = string
  description = "Image name"
}

variable task_role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
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

variable data_load_path {
  type        = string
  description = "S3 location to pull sql file from to load dev db data"
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
}
