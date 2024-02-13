data aws_region current {}

resource aws_ecs_task_definition task_definition {
  family        = "dp-${var.deployment_stage}-${var.custom_stack_name}-deletion"
  network_mode  = "awsvpc"
  cpu    = 2048
  memory = 4096
  task_role_arn = var.task_role_arn
  execution_role_arn = var.execution_role
  requires_compatibilities = [ "FARGATE" ]
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
      },
      {
        "name": "DOWNLOAD_WMG_DATA_TO_DISK",
        "value": "false"
      }
    ],
    "logConfiguration": {
      "logDriver": "awslogs",
      "options": {
        "awslogs-stream-prefix": "fargate",
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
