variable "deployment_stage" {
  type = string
}

// RDS

variable "db_username" {
  type = string
}
variable "db_password" {
  type = string
}
variable "db_instance_count" {
  type = string
}
variable "preferred_maintenance_window" {
  type = string
}

variable "log_bucket" {
  type = string
}
