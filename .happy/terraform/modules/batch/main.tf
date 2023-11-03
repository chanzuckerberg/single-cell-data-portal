# This is a batch job
# 

data aws_region current {}

data aws_caller_identity current {}

resource aws_batch_job_definition batch_job_def {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-upload"
  container_properties = jsonencode({
  "jobRoleArn": "${var.batch_role_arn}",
  "image": "${var.image}",
  "memory": var.batch_container_memory_limit,
  "environment": [
    {
      "name": "ARTIFACT_BUCKET",
      "value": "${var.artifact_bucket}"
    },
    {
      "name": "CELLXGENE_BUCKET",
      "value": "${var.cellxgene_bucket}"
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
      "name": "AWS_DEFAULT_REGION",
      "value": "${data.aws_region.current.name}"
    },
    {
      "name": "REMOTE_DEV_PREFIX",
      "value": "${var.remote_dev_prefix}"
    },
    {
      "name": "FRONTEND_URL",
      "value": "${var.frontend_url}"
    }
  ],
  "vcpus": 8,
  "linuxParameters": {
     "maxSwap": 800000,
     "swappiness": 60
  },
  "logConfiguration": {
    "logDriver": "awslogs",
    "options": {
      "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
      "awslogs-region": "${data.aws_region.current.name}"
    }
  }
})
}

resource aws_batch_job_definition dataset_metadata_update {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-dataset-metadata-update"
  container_properties = jsonencode({
  "command": [\"python3\",
    \"-m\",
    \"backend.layers.processing.dataset_metadata_update\",
  ],
  "jobRoleArn": "${var.batch_role_arn}",
  "image": "${var.image}",
  "memory": var.batch_container_memory_limit,
  "environment": [
    {
      "name": "ARTIFACT_BUCKET",
      "value": "${var.artifact_bucket}"
    },
    {
      "name": "CELLXGENE_BUCKET",
      "value": "${var.cellxgene_bucket}"
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
      "name": "AWS_DEFAULT_REGION",
      "value": "${data.aws_region.current.name}"
    },
    {
      "name": "REMOTE_DEV_PREFIX",
      "value": "${var.remote_dev_prefix}"
    }
  ],
  "vcpus": 8,
  "linuxParameters": {
     "maxSwap": 800000,
     "swappiness": 60
  },
  "logConfiguration": {
    "logDriver": "awslogs",
    "options": {
      "awslogs-group": "${aws_cloudwatch_log_group.cloud_watch_logs_group.id}",
      "awslogs-region": "${data.aws_region.current.name}"
    }
  }
})
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/upload"
}
