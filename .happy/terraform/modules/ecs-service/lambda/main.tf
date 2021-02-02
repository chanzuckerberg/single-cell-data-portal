# This is a lambda job
# 
resource aws_lambda_function lambda_job_def {
  role          = var.lambda_execution_role
  function_name = "${var.custom_stack_name}-lambda"
  package_type  = "Image"
  image_uri     = var.image
  environment {
    variables = {
      ARTIFACT_BUCKET = var.artifact_bucket,
      CELLXGENE_BUCKET = var.cellxgene_bucket,
      DEPLOYMENT_STAGE = var.deployment_stage,
      REMOTE_DEV_PREFIX = "/${var.custom_stack_name}",
    }
  }
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "${var.custom_stack_name}/lambda-handler-error"
}
