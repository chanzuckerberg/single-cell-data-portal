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

variable db_dump_s3_uri {
  type        = string
  description = "S3 location where a pg_dump file is located for dumping from one env and restoring to another"
  default     = ""
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
}
