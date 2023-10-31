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

variable max_attempts {
  type        = number
  description = "The max number of attempts to run the Job before failing out"
  default     = 3
}

# The following variables are for the batch job definition
variable image {
  type        = string
  description = "Image name"
}

variable batch_role_arn {
  type        = string
  description = "ARN for the role assumed by tasks"
}

variable remote_dev_prefix {
  type        = string
  description = "S3 storage path / db schema prefix"
  default     = ""
}

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

variable frontend_url {
  type        = string
  description = "url for the frontend app"
}

variable "batch_job_log_group" {
  type        = string
  description = "The name of the log group for the batch job"
}
