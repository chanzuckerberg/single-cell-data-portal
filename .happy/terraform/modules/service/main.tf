# This is a service managed by ECS attached to the environment's load balancer
#

data aws_region current {}

resource aws_ecs_service service {
  cluster         = var.cluster
  desired_count   = var.desired_count
  task_definition = aws_ecs_task_definition.task_definition.id
  launch_type     = "FARGATE"
  name            = "${var.custom_stack_name}-${var.app_name}"
  load_balancer {
    container_name   = "web"
    container_port   = var.service_port
    target_group_arn = aws_lb_target_group.target_group.id
  }
  network_configuration {
    security_groups  = var.security_groups
    subnets          = var.subnets
    assign_public_ip = false
  }
  deployment_circuit_breaker {
    enable   = true
    rollback = !var.fail_fast
  }
  enable_execute_command = true
  wait_for_steady_state = var.wait_for_steady_state
}

resource aws_ecs_task_definition task_definition {
  family                   = "dp-${var.deployment_stage}-${var.custom_stack_name}-${var.app_name}"
  memory                   = var.memory
  cpu                      = var.cpu
  network_mode             = "awsvpc"
  task_role_arn            = var.task_role_arn
  execution_role_arn       = var.execution_role
  requires_compatibilities = ["FARGATE"]
  container_definitions = <<EOF
[
  {
    "name": "web",
    "essential": true,
    "image": "${var.image}",
    "memory": ${var.memory},
    "environment": [
      {
        "name": "REMOTE_DEV_PREFIX",
        "value": "${var.remote_dev_prefix}"
      },
      {
        "name": "DATASET_SUBMISSIONS_BUCKET",
        "value": "${var.dataset_submissions_bucket}"
      },
      {
        "name": "DATASETS_BUCKET",
        "value": "${var.datasets_bucket}"
      },
      {
        "name": "DEPLOYMENT_STAGE",
        "value": "${var.deployment_stage}"
      },
      {
        "name": "AWS_REGION",
        "value": "${data.aws_region.current.name}"
      },
      {
        "name": "FRONTEND_URL",
        "value": "${var.frontend_url}"
      },
      {
        "name": "UPLOAD_SFN_ARN",
        "value": "${var.step_function_arn}"
      },
      {
        "name": "API_URL",
        "value": "${var.api_url}"
      },
      {
        "name": "AWS_DEFAULT_REGION",
        "value": "${data.aws_region.current.name}"
      }
    ],
    "portMappings": [
      {
        "containerPort": ${var.service_port}
      }
    ],
    "logConfiguration": {
      "logDriver": "awslogs",
      "options": {
        "awslogs-stream-prefix" : "fargate",
        "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
        "awslogs-region": "${data.aws_region.current.name}"
      }
    },
    "command": ${jsonencode((length(var.cmd) == 0) ? null : var.cmd)}
  }
]
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/${var.app_name}"
}

resource aws_lb_target_group target_group {
  vpc_id               = var.vpc
  port                 = var.service_port
  protocol             = "HTTP"
  target_type          = "ip"
  deregistration_delay = 10
  health_check {
    interval            = var.health_check_interval
    path                = "/"
    protocol            = "HTTP"
    timeout             = 5
    healthy_threshold   = 2
    unhealthy_threshold = 10
    matcher             = "200-299"
  }
}

resource aws_lb_listener_rule listener_rule {
  listener_arn = var.listener
  priority     = var.priority
  # Dev stacks need to match on hostnames
  dynamic "condition" {
    for_each = length(var.host_match) == 0 ? [] : [var.host_match]
    content {
      host_header {
        values = [
          condition.value
        ]
      }
    }
  }
  # Staging/prod envs are only expected to have a single stack,
  # so let's add all requests to that stack.
  dynamic "condition" {
    for_each = length(var.host_match) == 0 ? ["/*"] : []
    content {
      path_pattern {
        values = [condition.value]
      }
    }
  }
  action {
    target_group_arn = aws_lb_target_group.target_group.id
    type             = "forward"
  }
}
