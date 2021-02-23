resource "aws_sfn_state_machine" "state_machine" {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-sfn"
  role_arn = var.role_arn

  definition = <<EOF
{
    "StartAt": "Manage Batch task",
    "States": {
      "Manage Batch task": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Parameters": {
          "JobDefinition": "${var.job_definition_arn}",
          "JobName": "processing",
          "JobQueue": "${var.job_queue_arn}",
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "DROPBOX_URL",
                "Value.$": "$.url"
              },
              {
                "Name": "DATASET_ID",
                "Value.$": "$.dataset_uuid"
              }
            ]
          }
        },
        "End": true,
        "TimeoutSeconds": 10800,
        "Retry": [
          {
            "ErrorEquals": [
              "States.TaskFailed"
            ],
            "IntervalSeconds": 1,
            "BackoffRate": 2,
            "MaxAttempts": 2
          }
        ],
        "Catch": [
          {
            "ErrorEquals": [
              "States.ALL"
            ],
            "Next": "HandleErrors",
            "ResultPath": "$.error"
          }
        ]
      },
      "HandleErrors": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_error_handler}",
        "End": true
      }
    }
}
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/upload-sfn"
}
