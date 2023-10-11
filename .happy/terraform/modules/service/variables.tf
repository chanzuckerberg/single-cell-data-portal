variable "vpc" {
  type        = string
  description = "The VPC that the ECS cluster is deployed to"
}

variable "custom_stack_name" {
  type        = string
  description = "Please provide the stack name"
}

variable "remote_dev_prefix" {
  type        = string
  description = "S3 storage path / db schema prefix"
}

variable "dataset_submissions_bucket" {
  type        = string
  description = "S3 bucket for dataset submissions"
  default     = ""
}

variable datasets_bucket {
  type        = string
  description = "Datasets public-access bucket name"
  default     = ""
}

variable "app_name" {
  type        = string
  description = "Please provide the ECS service name"
}

variable "cluster" {
  type        = string
  description = "Please provide the ECS Cluster ID that this service should run on"
}

variable "image" {
  type        = string
  description = "Image name"
}

variable "service_port" {
  type        = number
  description = "What port does this service run on?"
  default     = 80
}

variable "desired_count" {
  type        = number
  description = "How many instances of this task should we run across our cluster?"
  default     = 2
}

variable "listener" {
  type        = string
  description = "The Application Load Balancer listener to register with"
}

variable "host_match" {
  type        = string
  description = "Host header to match for target rule. Leave empty to match all requests"
}

variable "host_match_internal" {
  type        = string
  description = "Host header to match for target rule. Leave empty to match all requests"
  default     = ""
}

variable "security_groups" {
  type        = list(string)
  description = "Security groups for ECS tasks"
}

variable "subnets" {
  type        = list(string)
  description = "Subnets for ecs tasks"
}

variable "task_role_arn" {
  type        = string
  description = "ARN for the role assumed by tasks"
}

variable "path" {
  type        = string
  description = "The path to register with the Application Load Balancer"
  default     = "/*"
}

variable "cmd" {
  type        = list(string)
  description = "The path to register with the Application Load Balancer"
  default     = []
}

variable "api_url" {
  type        = string
  description = "URL for the backend api."
}

variable "frontend_url" {
  type        = string
  description = "URL for the frontend app."
}

variable "api_domain" {
  type        = string
  description = "domain for the backend api"
}

variable "web_domain" {
  type        = string
  description = "domain for the website"
}

variable "data_locator_domain" {
  type        = string
  description = "domain for the data portal"
}

variable "cxg_bucket_path" {
  type        = string
  description = "path to the cxg bucket"
}

variable "deployment_stage" {
  type        = string
  description = "The name of the deployment stage of the Application"
  default     = "test"
}

variable "step_function_arn" {
  type        = string
  description = "ARN for the step function called by the uploader"
}

variable "priority" {
  type        = number
  description = "Listener rule priority number within the given listener"
}

variable "memory" {
  type        = number
  description = "Amount of memory to allocate to each task"
  default     = 2048
}

variable "cpu" {
  type        = number
  description = "CPU shares (1cpu=1024) per task"
  default     = 2048
}

variable "wait_for_steady_state" {
  type        = bool
  description = "Whether Terraform should block until the service is in a steady state before exiting"
  default     = false
}

variable "execution_role" {
  type        = string
  description = "Execution role to use for fargate tasks - required for fargate services!"
}

variable "health_check_interval" {
  type        = number
  description = "Interval for the health check pings"
  default     = 15
}

variable dd_key_secret_arn {
  type        = string
  description = "ARN for the Datadog secret key"
}
