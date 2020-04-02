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

// Browser RDS

variable "browser_db_username" {
  type = string
}
variable "browser_db_password" {
  type = string
}
variable "browser_db_instance_count" {
  type = string
}
variable "browser_preferred_maintenance_window" {
  type = string
}

// browser frontend
variable "refer_secret" {
  type = string
}
