variable image {
  type        = string
  description = "Image name"
}

variable artifact_bucket {
  type        = string
  description = "Artifact bucket name"
}

variable cellxgene_bucket {
  type        = string
  description = "Cellxgene bucket name"
}

variable custom_stack_name {
  type        = string
  description = "Please provide the stack name"
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
}

variable lambda_execution_role {
  type        = string
  description = "Role for lambda execution"
}
