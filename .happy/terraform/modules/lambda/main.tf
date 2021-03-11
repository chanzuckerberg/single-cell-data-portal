# This is a lambda job
# 
resource aws_lambda_function lambda_job_def {
  role          = var.lambda_execution_role
  function_name = "dp-${var.deployment_stage}-${var.custom_stack_name}-${var.name}"
  package_type  = "Image"
  image_uri     = var.image
  timeout       = 60
  environment {
    variables = {
      ARTIFACT_BUCKET = var.artifact_bucket,
      CELLXGENE_BUCKET = var.cellxgene_bucket,
      DEPLOYMENT_STAGE = var.deployment_stage,
      REMOTE_DEV_PREFIX = var.remote_dev_prefix
    }
  }
  vpc_config {
    subnet_ids         = var.subnets
    security_group_ids = var.security_groups
  }
}
