variable job_queue_arn {
  type        = string
  description = "ARN of the batch job queue"
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

variable artifact_bucket {
  type        = string
  description = "Artifact bucket name"
}

variable datasets_bucket {
  type        = string
  description = "Datasets public-access bucket name"
}

variable image {
  type        = string
  description = "Image name"
}

variable batch_role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
}

variable cmd {
  type        = string
  description = "Command to run"
  default     = ""
}

variable remote_dev_prefix {
  type        = string
  description = "S3 storage path / db schema prefix"
  default     = ""
}

variable batch_container_memory_limit {
  type        = number
  description = "Memory hard limit for the batch container"
  default     = 28000
}

variable sfn_role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
}