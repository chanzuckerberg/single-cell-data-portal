variable artifact_bucket {
  type        = string
  description = "Artifact bucket name"
}

variable wmg_bucket {
  type        = string
  description = "Where's my gene (output) bucket name"
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

variable batch_container_memory_limit {
  type        = number
  description = "Memory hard limit for the batch container"
  default     = 28000
}

variable desired_vcpus {
  type        = number
  description = "Number of desired vCPUs"
  default     = 2
}

variable "api_url" {
  type        = string
  description = "URL for the backend api."
}