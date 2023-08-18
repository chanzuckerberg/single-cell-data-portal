data aws_region current {}

data aws_caller_identity current {}

resource aws_batch_job_definition batch_job_def {
  type = "container"
  name = "wmg-${var.deployment_stage}-${var.custom_stack_name}-wmg-processing"
  container_properties = <<EOF
{
  "jobRoleArn": "${var.batch_role_arn}",
  "image": "${var.image}",
  "memory": ${var.batch_container_memory_limit},
  "environment": [
    {
      "name": "ARTIFACT_BUCKET",
      "value": "${var.artifact_bucket}"
    },
    {
      "name": "WMG_BUCKET",
      "value": "${var.wmg_bucket}"
    },
    {
      "name": "API_URL",
      "value": "${var.api_url}"
    },
    {
      "name": "DEPLOYMENT_STAGE",
      "value": "${var.deployment_stage}"
    },
    {
      "name": "AWS_DEFAULT_REGION",
      "value": "${data.aws_region.current.name}"
    },
    {
      "name": "REMOTE_DEV_PREFIX",
      "value": "${var.remote_dev_prefix}"
    }
  ],
  "vcpus": ${var.desired_vcpus},
  "logConfiguration": {
    "logDriver": "awslogs",
    "options": {
      "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
      "awslogs-region": "${data.aws_region.current.name}"
    }
  }
}
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/wmg-processing"
}