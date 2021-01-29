# This deploys a dev stack.
# 
module frontend_dns {
  source                = "../dns"
  custom_stack_name     = var.custom_stack_name
  app_name              = "frontend"
  alb_dns               = var.frontend_alb_dns
  canonical_hosted_zone = var.frontend_alb_zone
  zone                  = var.zone
}

module backend_dns {
  source                = "../dns"
  custom_stack_name     = var.custom_stack_name
  app_name              = "backend"
  alb_dns               = var.backend_alb_dns
  canonical_hosted_zone = var.backend_alb_zone
  zone                  = var.zone
}

module frontend_service {
  source            = "../service"
  custom_stack_name = var.custom_stack_name
  app_name          = "frontend"
  vpc               = var.vpc
  image             = "${var.frontend_image_repo}:${var.image_tag}"
  cluster           = var.cluster
  desired_count     = 2
  listener          = var.frontend_listener_arn
  subnets           = var.subnets
  security_groups   = var.security_groups
  task_role_arn     = var.task_role_arn
  service_port      = 9000
  deployment_stage  = var.deployment_stage
  step_function_arn = module.upload_sfn.step_function_arn
  host_match        = join(".", [module.frontend_dns.dns_prefix, var.external_dns])
  priority          = var.priority
  api_url           = join("", ["https://", module.backend_dns.dns_prefix, ".", var.external_dns])
  frontend_url      = join("", ["https://", module.frontend_dns.dns_prefix, ".", var.external_dns])
}

module backend_service {
  source            = "../service"
  custom_stack_name = var.custom_stack_name
  app_name          = "backend"
  vpc               = var.vpc
  image             = "${var.backend_image_repo}:${var.image_tag}"
  cluster           = var.cluster
  desired_count     = 2
  listener          = var.backend_listener_arn
  subnets           = var.subnets
  security_groups   = var.security_groups
  task_role_arn     = var.task_role_arn
  service_port      = 5000
  cmd               = var.backend_cmd
  deployment_stage  = var.deployment_stage
  step_function_arn = module.upload_sfn.step_function_arn
  host_match        = join(".", [module.backend_dns.dns_prefix, var.external_dns])
  priority          = var.priority
  api_url           = join("", ["https://", module.backend_dns.dns_prefix, ".", var.external_dns])
  frontend_url      = join("", ["https://", module.frontend_dns.dns_prefix, ".", var.external_dns])
}

module migrate_db {
  source            = "../migration"
  image             = "${var.backend_image_repo}:${var.image_tag}"
  task_role_arn     = var.task_role_arn
  cmd               = var.migration_cmd
  custom_stack_name = var.custom_stack_name
  deployment_stage  = var.deployment_stage
  data_load_path    = var.data_load_path
}

module delete_db {
  source            = "../deletion"
  image             = "${var.backend_image_repo}:${var.image_tag}"
  task_role_arn     = var.task_role_arn
  cmd               = var.deletion_cmd
  custom_stack_name = var.custom_stack_name
  deployment_stage  = var.deployment_stage
}

module upload_batch {
  source            = "../batch"
  image             = "${var.upload_image_repo}:${var.image_tag}"
  batch_role_arn    = var.batch_role_arn
  cmd               = ""
  custom_stack_name = var.custom_stack_name
  deployment_stage  = var.deployment_stage
  artifact_bucket   = var.artifact_bucket
  cellxgene_bucket  = var.cellxgene_bucket
  frontend_url      = join("", ["https://", module.frontend_dns.dns_prefix, ".", var.external_dns])
}

module upload_lambda {
  source                = "../lambda"
  image                 = "${var.lambda_upload_repo}:${var.image_tag}"
  custom_stack_name     = var.custom_stack_name
  deployment_stage      = var.deployment_stage
  artifact_bucket       = var.artifact_bucket
  cellxgene_bucket      = var.cellxgene_bucket
  lambda_execution_role = var.lambda_execution_role
}

module upload_sfn {
  source               = "../sfn"
  job_definition_arn   = module.upload_batch.batch_job_definition
  job_queue_arn        = var.job_queue_arn
  role_arn             = var.sfn_role_arn
  custom_stack_name    = var.custom_stack_name
  lambda_error_handler = module.upload_lambda.error_handler
}
