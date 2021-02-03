variable custom_stack_name {
  type        = string
  description = "Please provide the stack name"
}

variable vpc {
  type        = string
  description = "Choose which VPC the Application Load Balancer should be deployed to"
}

variable cluster {
  type        = string
  description = "Which cluster to launch ECS tasks in"
}

variable subnets {
  type        = string
  description = "Choose which subnets the Application Load Balancer should be deployed to"
}

variable security_groups {
  type        = string
  description = "Select the Security Groups to apply to the Application Load Balancer"
}

variable task_role_arn {
  type        = string
  description = "Role ARN for ECS tasks"
}

variable batch_role_arn {
  type        = string
  description = "Role ARN for Batch Jobs"
}

variable backend_listener_arn {
  type        = string
  description = "ARN for backend alb listener"
}

variable frontend_listener_arn {
  type        = string
  description = "ARN for frontend alb listener"
}

variable backend_alb_dns {
  type        = string
  description = "DNS name for backend alb"
}

variable frontend_alb_dns {
  type        = string
  description = "DNS name for frontend alb"
}

variable backend_alb_zone {
  type        = string
  description = "R53 zone for backend alb"
}

variable frontend_alb_zone {
  type        = string
  description = "R53 zone for frontend alb"
}

variable backend_image_repo {
  type        = string
  description = "Docker image to use for the backend"
}

variable frontend_image_repo {
  type        = string
  description = "Docker image to use for the frontend"
}

variable image_tag {
  type        = string
  description = "Tagged image to use for this remote-dev env."
}

variable zone {
  type        = string
  description = "Route53 zone name"
}

variable external_dns {
  type        = string
  description = "Domain suffix for our multidomain oauth proxy"
}

variable backend_cmd {
  type        = string
  description = "Entrypoint cmd for frontend containers"
}

variable frontend_cmd {
  type        = string
  description = "Entrypoint cmd for frontend containers"
}

variable migration_cmd {
  type        = string
  description = "Command to run when this ECS service is updated."
}

variable deletion_cmd {
  type        = string
  description = "Command to run when this ECS service is updated."
}

variable deployment_stage {
  type        = string
  description = "The name of the deployment stage of the Application"
}

variable priority {
  type        = number
  description = "Listener rule priority number within the given listener"
}

variable job_queue_arn {
  type        = string
  description = "ARN of the batch job queue"
}

variable upload_image_repo {
  type        = string
  description = "Docker image to use for the upload batch job"
}

variable sfn_role_arn {
  type        = string
  description = "Role ARN for Step Functions"
}

variable artifact_bucket {
  type        = string
  description = "S3 bucket for upload artifacts"
}

variable cellxgene_bucket {
  type        = string
  description = "S3 bucket for cellxgene outputs"
}

variable data_load_path {
  type        = string
  description = "S3 location to pull sql file from to load dev db data"
}

variable lambda_execution_role {
  type        = string
  description = "Role for lambda execution"
}

variable lambda_upload_repo {
  type        = string
  description = "Docker image to use for the lambda upload job"
}
