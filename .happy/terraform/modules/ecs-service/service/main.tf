# This is a service managed by ECS with an internally-facing load balancer. Public access to the load balancer needs to be managed by the multioauth proxy.
# 

data aws_region current {}

resource aws_ecs_service service {
  cluster         = var.cluster
  desired_count   = var.desired_count
  task_definition = aws_ecs_task_definition.task_definition.id
  launch_type     = "EC2"
  name            = "${var.custom_stack_name}-${var.app_name}"
  load_balancer {
    container_name   = "web"
    container_port   = var.service_port
    target_group_arn = aws_lb_target_group.target_group.id
  }
  network_configuration {
    security_groups  = split(",", var.security_groups)
    subnets          = split(",", var.subnets)
    assign_public_ip = false
  }
}

resource aws_ecs_task_definition task_definition {
  family        = "${var.custom_stack_name}-${var.app_name}"
  network_mode  = "awsvpc"
  task_role_arn = var.task_role_arn
  container_definitions = <<EOF
[
  {
    "name": "web",
    "essential": true,
    "image": "${var.image}",
    "memory": 512,
    "environment": [
      {
        "name": "REMOTE_DEV_PREFIX",
        "value": "/${var.custom_stack_name}"
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
        "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
        "awslogs-region": "${data.aws_region.current.name}"
      }
    },
    "command": ${jsonencode(!(var.cmd == "") ? split(",", var.cmd) : null)}
  }
]
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "${var.custom_stack_name}/${var.app_name}"
}

resource aws_lb_target_group target_group {
  vpc_id               = var.vpc
  port                 = var.service_port
  protocol             = "HTTP"
  target_type          = "ip"
  deregistration_delay = 10
  health_check {
    interval            = 15
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
  condition {
    host_header {
      values = [
        var.host_match,
      ]
    }
  }
  action {
    target_group_arn = aws_lb_target_group.target_group.id
    type             = "forward"
  }
}


