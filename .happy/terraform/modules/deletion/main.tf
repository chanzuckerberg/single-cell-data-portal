data aws_region current {}

resource aws_ecs_task_definition task_definition {
  family        = "dp-${var.deployment_stage}-${var.custom_stack_name}-deletion"
  network_mode  = "awsvpc"
  task_role_arn = var.task_role_arn
  container_definitions = <<EOF
[
  {
    "name": "deletedb",
    "essential": true,
    "image": "${var.image}",
    "memory": 512,
    "environment": [
      {
        "name": "AWS_REGION",
        "value": "${data.aws_region.current.name}"
      },
      {
        "name": "AWS_DEFAULT_REGION",
        "value": "${data.aws_region.current.name}"
      },
      {
        "name": "REMOTE_DEV_PREFIX",
        "value": "${var.remote_dev_prefix}"
      },
      {
        "name": "DEPLOYMENT_STAGE",
        "value": "${var.deployment_stage}"
      }
    ],
    "logConfiguration": {
      "logDriver": "awslogs",
      "options": {
        "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
        "awslogs-region": "${data.aws_region.current.name}"
      }
    },
    "command": ${jsonencode(var.cmd)}
  }
]
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/deletion"
}
