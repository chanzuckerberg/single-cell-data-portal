variable artifact_bucket {
  type        = string
  description = "Artifact bucket name"
}

variable cellxgene_bucket {
  type        = string
  description = "Cellxgene bucket name"
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

variable frontend_url {
  type        = string
  description = "url for the frontend app"
}

variable batch_container_memory_limit {
  type        = number
  description = "Memory hard limit for the batch container"
  default     = 63500
}
