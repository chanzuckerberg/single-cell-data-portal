# This deploys a Data Portal stack.
# 

data aws_secretsmanager_secret_version config {
  secret_id = var.happy_config_secret
}

locals {
  secret = jsondecode(nonsensitive(data.aws_secretsmanager_secret_version.config.secret_string))
  alb_key = var.require_okta ? "private_albs" : "public_albs"

  custom_stack_name            = var.stack_name
  image_tag                    = var.image_tag
  priority                     = var.priority
  deployment_stage             = var.deployment_stage
  remote_dev_prefix            = var.stack_prefix
  wait_for_steady_state        = var.wait_for_steady_state
  batch_container_memory_limit = var.batch_container_memory_limit
  wmg_batch_container_memory_limit = var.wmg_batch_container_memory_limit
  wmg_desired_vcpus                = var.wmg_desired_vcpus
  cg_desired_vcpus                = var.cg_desired_vcpus
  cg_batch_container_memory_limit = var.cg_batch_container_memory_limit

  migration_cmd                = ["make", "-C", "/single-cell-data-portal/backend", "db/init_remote_dev"]
  deletion_cmd                 = ["make", "-C", "/single-cell-data-portal/backend", "db/delete_remote_dev"]
  frontend_cmd                 = []
  # TODO: Assess whether this is safe for Portal API as well. Trying 1 worker in rdev portal backend containers, to minimize use of memory by TileDB (allocates multi-GB per process)
  # Note: keep-alive timeout should always be greater than the idle timeout of the load balancer (60 seconds)
  backend_workers              = var.backend_workers 
  backend_cmd                  = ["gunicorn", "--worker-class", "gevent", "--workers", "${local.backend_workers}",
    "--bind", "0.0.0.0:5000", "backend.api_server.app:app", "--max-requests", "10000", "--timeout", "180",
    "--keep-alive", "61", "--log-level", "info"]
  data_load_path               = "s3://${local.secret["s3_buckets"]["env"]["name"]}/database/seed_data_04_2f30f3bcc9aa.sql"

  vpc_id                          = local.secret["cloud_env"]["vpc_id"]
  subnets                         = local.secret["cloud_env"]["private_subnets"]
  security_groups                 = local.secret["security_groups"]
  zone                            = local.secret["zone_id"]
  cluster                         = local.secret["cluster_arn"]
  ecs_execution_role              = lookup(local.secret, "ecs_execution_role", "")

  frontend_image_repo             = local.secret["ecrs"]["frontend"]["url"]
  backend_image_repo              = local.secret["ecrs"]["backend"]["url"]
  upload_image_repo               = local.secret["ecrs"]["processing"]["url"]
  lambda_upload_success_repo      = local.secret["ecrs"]["upload_success"]["url"]
  lambda_upload_repo              = local.secret["ecrs"]["upload_failures"]["url"]
  lambda_dataset_submissions_repo = local.secret["ecrs"]["dataset_submissions"]["url"]
  wmg_upload_image_repo           = local.secret["ecrs"]["wmg_processing"]["url"]
  cg_upload_image_repo            = local.secret["ecrs"]["cellguide_pipeline"]["url"]
  batch_role_arn                  = local.secret["batch_queues"]["upload"]["role_arn"]
  job_queue_arn                   = local.secret["batch_queues"]["upload"]["queue_arn"]
  wmg_batch_role_arn              = local.secret["batch_queues"]["wmg"]["role_arn"]
  cg_batch_role_arn               = local.secret["batch_queues"]["cg"]["role_arn"]
  schema_migration_job_queue_arn  = local.secret["batch_queues"]["schema_migration"]["queue_arn"]
  schema_migration_batch_role_arn = local.secret["batch_queues"]["schema_migration"]["role_arn"]
  external_dns                    = local.secret["external_zone_name"]
  internal_dns                    = local.secret["internal_zone_name"]

  frontend_listener_arn        = try(local.secret[local.alb_key]["frontend"]["listener_arn"], "")
  backend_listener_arn         = try(local.secret[local.alb_key]["backend"]["listener_arn"], "")
  frontend_alb_zone            = try(local.secret[local.alb_key]["frontend"]["zone_id"], "")
  backend_alb_zone             = try(local.secret[local.alb_key]["backend"]["zone_id"], "")
  frontend_alb_dns             = try(local.secret[local.alb_key]["frontend"]["dns_name"], "")
  backend_alb_dns              = try(local.secret[local.alb_key]["backend"]["dns_name"], "")

  artifact_bucket            = try(local.secret["s3_buckets"]["artifact"]["name"], "")
  cellxgene_bucket           = try(local.secret["s3_buckets"]["cellxgene"]["name"], "")
  datasets_bucket            = try(local.secret["s3_buckets"]["datasets"]["name"], "")
  dataset_submissions_bucket = try(local.secret["s3_buckets"]["dataset_submissions"]["name"], "")
  wmg_bucket                 = try(local.secret["s3_buckets"]["wmg"]["name"], "")

  ecs_role_arn                 = local.secret["service_roles"]["ecs_role"]
  sfn_role_arn                 = local.secret["service_roles"]["sfn_upload"]
  lambda_execution_role        = local.secret["service_roles"]["lambda_errorhandler"]
  schema_migration_sfn_role_arn = local.secret["service_roles"]["schema_migration_sfn_service"]

  frontend_url = try(join("", ["https://", module.frontend_dns[0].dns_prefix, ".", local.external_dns]), var.frontend_url)
  backend_url  = try(join("", ["https://", module.backend_dns[0].dns_prefix, ".", local.external_dns]), var.backend_url)
}

module frontend_dns {
  count                 = var.require_okta ? 1 : 0
  source                = "../dns"
  custom_stack_name     = local.custom_stack_name
  app_name              = "frontend"
  alb_dns               = local.frontend_alb_dns
  canonical_hosted_zone = local.frontend_alb_zone
  zone                  = local.internal_dns
}

module backend_dns {
  count                 = var.require_okta ? 1 : 0
  source                = "../dns"
  custom_stack_name     = local.custom_stack_name
  app_name              = "backend"
  alb_dns               = local.backend_alb_dns
  canonical_hosted_zone = local.backend_alb_zone
  zone                  = local.internal_dns
}

module frontend_service {
  source                     = "../service"
  custom_stack_name          = local.custom_stack_name
  app_name                   = "frontend"
  vpc                        = local.vpc_id
  image                      = "${local.frontend_image_repo}:${local.image_tag}"
  cluster                    = local.cluster
  desired_count              = var.frontend_instance_count
  listener                   = local.frontend_listener_arn
  subnets                    = local.subnets
  security_groups            = local.security_groups
  task_role_arn              = local.ecs_role_arn
  service_port               = 9000
  memory                     = var.frontend_memory
  cpu                        = 2048
  deployment_stage           = local.deployment_stage
  step_function_arn          = module.upload_sfn.step_function_arn
  host_match                 = try(join(".", [module.frontend_dns[0].dns_prefix, local.external_dns]), "")
  priority                   = local.priority
  api_url                    = local.backend_url
  frontend_url               = local.frontend_url
  remote_dev_prefix          = local.remote_dev_prefix
  dataset_submissions_bucket = local.dataset_submissions_bucket
  execution_role             = local.ecs_execution_role
  health_check_interval      = 15
  dd_key_secret_arn          = var.dd_key_secret_arn

  wait_for_steady_state = local.wait_for_steady_state
}

module backend_service {
  source                     = "../service"
  custom_stack_name          = local.custom_stack_name
  app_name                   = "backend"
  vpc                        = local.vpc_id
  image                      = "${local.backend_image_repo}:${local.image_tag}"
  cluster                    = local.cluster
  desired_count              = var.backend_instance_count
  listener                   = local.backend_listener_arn
  subnets                    = local.subnets
  security_groups            = local.security_groups
  task_role_arn              = local.ecs_role_arn
  service_port               = 5000
  memory                     = var.backend_memory
  cpu                        = var.backend_cpus * 1024
  cmd                        = local.backend_cmd
  deployment_stage           = local.deployment_stage
  step_function_arn          = module.upload_sfn.step_function_arn
  host_match                 = try(join(".", [module.backend_dns[0].dns_prefix, local.external_dns]), "")
  priority                   = local.priority
  api_url                    = local.backend_url
  frontend_url               = local.frontend_url
  remote_dev_prefix          = local.remote_dev_prefix
  dataset_submissions_bucket = local.dataset_submissions_bucket
  datasets_bucket            = local.datasets_bucket
  execution_role             = local.ecs_execution_role
  dd_key_secret_arn          = var.dd_key_secret_arn

  wait_for_steady_state = local.wait_for_steady_state
}

module migrate_db {
  source            = "../migration"
  image             = "${local.backend_image_repo}:${local.image_tag}"
  task_role_arn     = local.ecs_role_arn
  cmd               = local.migration_cmd
  custom_stack_name = local.custom_stack_name
  remote_dev_prefix = local.remote_dev_prefix
  deployment_stage  = local.deployment_stage
  data_load_path    = local.data_load_path
  execution_role    = local.ecs_execution_role
}

module delete_db {
  count             = var.delete_protected ? 0 : 1
  source            = "../deletion"
  image             = "${local.backend_image_repo}:${local.image_tag}"
  task_role_arn     = local.ecs_role_arn
  cmd               = local.deletion_cmd
  custom_stack_name = local.custom_stack_name
  remote_dev_prefix = local.remote_dev_prefix
  deployment_stage  = local.deployment_stage
  execution_role    = local.ecs_execution_role
}

module upload_batch {
  source            = "../batch"
  image             = "${local.upload_image_repo}:${local.image_tag}"
  batch_role_arn    = local.batch_role_arn
  cmd               = ""
  custom_stack_name = local.custom_stack_name
  remote_dev_prefix = local.remote_dev_prefix
  deployment_stage  = local.deployment_stage
  artifact_bucket   = local.artifact_bucket
  cellxgene_bucket  = local.cellxgene_bucket
  datasets_bucket   = local.datasets_bucket
  frontend_url      = local.frontend_url
  batch_container_memory_limit = local.batch_container_memory_limit
}

module wmg_batch {
  source                        = "../wmg-batch"
  image                         = "${local.wmg_upload_image_repo}:${local.image_tag}"
  batch_role_arn                = local.wmg_batch_role_arn
  cmd                           = ""
  api_url                       = local.backend_url
  custom_stack_name             = local.custom_stack_name
  remote_dev_prefix             = local.remote_dev_prefix
  deployment_stage              = local.deployment_stage
  artifact_bucket               = local.artifact_bucket
  wmg_bucket                    = local.wmg_bucket
  desired_vcpus                 = local.wmg_desired_vcpus
  batch_container_memory_limit  = var.wmg_batch_container_memory_limit
}

module cg_batch {
  source                        = "../cg-batch"
  image                         = "${local.cg_upload_image_repo}:${local.image_tag}"
  batch_role_arn                = local.cg_batch_role_arn
  cmd                           = ""
  api_url                       = local.backend_url
  custom_stack_name             = local.custom_stack_name
  remote_dev_prefix             = local.remote_dev_prefix
  deployment_stage              = local.deployment_stage
  artifact_bucket               = local.artifact_bucket
  wmg_bucket                    = local.wmg_bucket
  desired_vcpus                 = local.cg_desired_vcpus
  batch_container_memory_limit  = var.cg_batch_container_memory_limit
}


module upload_success_lambda {
  source                     = "../lambda"
  image                      = "${local.lambda_upload_success_repo}:${local.image_tag}"
  name                       = "upload-success"
  custom_stack_name          = local.custom_stack_name
  remote_dev_prefix          = local.remote_dev_prefix
  deployment_stage           = local.deployment_stage
  artifact_bucket            = local.artifact_bucket
  cellxgene_bucket           = local.cellxgene_bucket
  lambda_execution_role      = local.lambda_execution_role
  subnets                    = local.subnets
  security_groups            = local.security_groups
  datasets_bucket            = local.datasets_bucket
}

module upload_error_lambda {
  source                     = "../lambda"
  image                      = "${local.lambda_upload_repo}:${local.image_tag}"
  name                       = "uploadfailures"
  custom_stack_name          = local.custom_stack_name
  remote_dev_prefix          = local.remote_dev_prefix
  deployment_stage           = local.deployment_stage
  artifact_bucket            = local.artifact_bucket
  cellxgene_bucket           = local.cellxgene_bucket
  datasets_bucket            = local.datasets_bucket
  lambda_execution_role      = local.lambda_execution_role
  subnets                    = local.subnets
  security_groups            = local.security_groups
}

module upload_sfn {
  source                 = "../sfn"
  job_definition_arn     = module.upload_batch.batch_job_definition_no_revision
  job_queue_arn          = local.job_queue_arn
  role_arn               = local.sfn_role_arn
  custom_stack_name      = local.custom_stack_name
  lambda_success_handler = module.upload_success_lambda.arn
  lambda_error_handler   = module.upload_error_lambda.arn
  deployment_stage       = local.deployment_stage

  # use for defining the batch job definition
  image             = "${local.upload_image_repo}:${local.image_tag}"
  batch_role_arn    = local.batch_role_arn
  remote_dev_prefix = local.remote_dev_prefix
  artifact_bucket   = local.artifact_bucket
  cellxgene_bucket  = local.cellxgene_bucket
  datasets_bucket   = local.datasets_bucket
  frontend_url      = local.frontend_url
  batch_job_log_group = module.upload_batch.batch_job_log_group
}

module dataset_submissions_lambda {
  source                     = "../lambda"
  image                      = "${local.lambda_dataset_submissions_repo}:${local.image_tag}"
  name                       = "dataset-submissions"
  custom_stack_name          = local.custom_stack_name
  remote_dev_prefix          = local.remote_dev_prefix
  deployment_stage           = local.deployment_stage
  artifact_bucket            = local.artifact_bucket
  cellxgene_bucket           = local.cellxgene_bucket
  datasets_bucket            = local.datasets_bucket
  lambda_execution_role      = aws_iam_role.dataset_submissions_lambda_service_role.arn
  step_function_arn          = module.upload_sfn.step_function_arn
  subnets                    = local.subnets
  security_groups            = local.security_groups
}

module schema_migration {
  source                        = "../schema_migration"
  image                         = "${local.upload_image_repo}:${local.image_tag}"
  batch_role_arn                = local.schema_migration_batch_role_arn
  cmd                           = ""
  custom_stack_name             = local.custom_stack_name
  remote_dev_prefix             = local.remote_dev_prefix
  deployment_stage              = local.deployment_stage
  artifact_bucket               = local.artifact_bucket
  datasets_bucket               = local.datasets_bucket
  batch_container_memory_limit  = local.batch_container_memory_limit
  job_queue_arn                 = local.schema_migration_job_queue_arn
  sfn_role_arn                  = local.schema_migration_sfn_role_arn

}

resource "aws_iam_role" "dataset_submissions_lambda_service_role" {
  name               = "corpora-dataset-submissions-service-role-${local.custom_stack_name}"
  path               = "/service-role/"
  assume_role_policy = data.aws_iam_policy_document.assume_role.json
}

data "aws_iam_policy_document" "assume_role" {
  statement {
    principals {
      type        = "Service"
      identifiers = ["lambda.amazonaws.com"]
    }

    actions = ["sts:AssumeRole"]
  }
}

data "aws_iam_policy_document" "lambda_step_function_execution_policy_document" {
  statement {
    sid    = "sfn"
    effect = "Allow"
    actions = [
      "states:StartExecution"
    ]
    resources = [
      module.upload_sfn.step_function_arn
    ]
  }
  statement {
    sid      = "cw"
    effect = "Allow"
    actions = [
      "logs:CreateLogGroup",
      "logs:CreateLogStream",
      "logs:PutLogEvents"
    ]
    resources = [
      "arn:aws:logs:*:*:*"
    ]
  }
  statement {
    sid    = "networking"
    effect = "Allow"
    actions = [
      "ec2:DescribeNetworkInterfaces",
      "ec2:CreateNetworkInterface",
      "ec2:DeleteNetworkInterface",
      "ec2:DescribeInstances",
      "ec2:AttachNetworkInterface",
      "ec2:DescribeSecurityGroups",
      "ec2:DescribeSubnets",
      "ec2:DescribeVpcs"
    ]
    resources = ["*"]
  }
  statement {
    actions = [
      "secretsmanager:DescribeSecret",
      "secretsmanager:GetSecretValue"
    ]
    resources = [
      "arn:aws:secretsmanager:us-west-2:${var.aws_account_id}:secret:corpora/backend/${var.deployment_stage}/*"
    ]
  }
  statement {
    sid     = "s3"
    effect  = "Allow"
    actions = [
      "s3:ListBucket",
      "s3:DeleteObject"
    ]
    resources = [
      "arn:aws:s3:::${local.artifact_bucket}",
      "arn:aws:s3:::${local.artifact_bucket}/*",
      "arn:aws:s3:::${local.cellxgene_bucket}",
      "arn:aws:s3:::${local.cellxgene_bucket}/*",
      "arn:aws:s3:::${local.datasets_bucket}",
      "arn:aws:s3:::${local.datasets_bucket}/*"
    ]
  }
}

resource "aws_iam_policy" "lambda_step_function_execution_policy" {
  name = "lambda-step-function-execution-policy-${local.custom_stack_name}"
  policy = data.aws_iam_policy_document.lambda_step_function_execution_policy_document.json
}

resource "aws_iam_role_policy_attachment" "lambda_step_function_execution_policy_attachment" {
  role       = aws_iam_role.dataset_submissions_lambda_service_role.name
  policy_arn = aws_iam_policy.lambda_step_function_execution_policy.arn
}

resource "aws_lambda_permission" "allow_dataset_submissions_lambda_execution" {
  statement_id  = "AllowExecutionFromS3Bucket"
  action        = "lambda:InvokeFunction"
  function_name = module.dataset_submissions_lambda.arn
  principal     = "s3.amazonaws.com"
  source_arn    = try(local.secret["s3_buckets"]["dataset_submissions"]["arn"], "")
}

resource "aws_s3_bucket_notification" "on_dataset_submissions_object_created" {
  bucket = local.dataset_submissions_bucket

  lambda_function {
    lambda_function_arn = module.dataset_submissions_lambda.arn
    events = ["s3:ObjectCreated:*"]
  }

  depends_on = [aws_lambda_permission.allow_dataset_submissions_lambda_execution]
}
