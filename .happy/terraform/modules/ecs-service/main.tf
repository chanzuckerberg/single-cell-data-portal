data aws_secretsmanager_secret_version config {
  secret_id = var.happy_config_secret
}

locals {
  secret = jsondecode(data.aws_secretsmanager_secret_version.config.secret_string)
}

module dev_env {
  source                = "./master"
  custom_stack_name     = var.stack_name
  vpc                   = local.secret["vpc"]
  subnets               = local.secret["subnets"]
  security_groups       = local.secret["security_groups"]
  zone                  = local.secret["zone"]
  external_dns          = local.secret["external_dns"]
  frontend_listener_arn = local.secret["frontend_listener_arn"]
  backend_listener_arn  = local.secret["backend_listener_arn"]
  frontend_alb_zone     = local.secret["frontend_alb_zone"]
  backend_alb_zone      = local.secret["backend_alb_zone"]
  frontend_alb_dns      = local.secret["frontend_alb_dns"]
  backend_alb_dns       = local.secret["backend_alb_dns"]
  frontend_image_repo   = local.secret["frontend_image_repo"]
  backend_image_repo    = local.secret["backend_image_repo"]
  image_tag             = var.image_tag
  task_role_arn         = local.secret["task_role_arn"]
  cluster               = local.secret["cluster_arn"]
  data_load_path        = "s3://${local.secret["env_s3_bucket"]}/database/dev_data.sql"
  migration_cmd         = "make,-C,/corpora-data-portal/backend,db/init_remote_dev"
  deletion_cmd          = "make,-C,/corpora-data-portal/backend,db/delete_remote_dev"
  frontend_cmd          = ""
  backend_cmd           = "python3,/chalice/run_local_server.py,--host,0.0.0.0"
  deployment_stage      = "rdev"
  priority              = var.priority
  sfn_role_arn          = local.secret["sfn_role_arn"]
  batch_role_arn        = local.secret["batch_role_arn"]
  upload_image_repo     = local.secret["upload_image_repo"]
  job_queue_arn         = local.secret["job_queue_arn"]
  artifact_bucket       = local.secret["artifact_bucket"]
  cellxgene_bucket      = local.secret["cellxgene_bucket"]
  lambda_upload_repo    = local.secret["lambda_upload_repo"]
  lambda_execution_role = local.secret["lambda_execution_role"]
}
