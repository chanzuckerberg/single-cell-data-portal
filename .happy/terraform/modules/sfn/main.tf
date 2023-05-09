# Same file as https://github.com/chanzuckerberg/single-cell-infra/blob/main/.happy/terraform/modules/sfn/main.tf
# This is used for environment (dev, staging, prod) deployments

resource "aws_sfn_state_machine" "state_machine" {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-sfn"
  role_arn = var.role_arn

  definition = <<EOF
{
    "StartAt": "DownloadValidate",
    "States": {
      "DownloadValidate": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Next": "CxgSeuratParallel",
        "Parameters": {
          "JobDefinition": "${var.job_definition_arn}",
          "JobName": "download-validate",
          "JobQueue": "${var.job_queue_arn}",
          "RetryStrategy": {
            "Attempts": ${var.max_attempts},
            "EvaluateOnExit": [
              {
                "Action": "EXIT",
                "OnExitCode": "1"
              },
              {
                "Action": "RETRY",
                "OnExitCode": "*"
              }
            ]
          },
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "DROPBOX_URL",
                "Value.$": "$.url"
              },
              {
                "Name": "DATASET_ID",
                 "Value.$": "$.dataset_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "download-validate"
              }
            ]
          }
        },
        "ResultPath": null,
        "TimeoutSeconds": 36000,
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
      "CxgSeuratParallel": {
        "Type": "Parallel",
        "Next": "HandleSuccess",
        "ResultPath": null,
        "Branches": [
          {
            "StartAt": "Cxg",
            "States": {
              "Cxg": {
                "Type": "Task",
                "End": true,
                "Resource": "arn:aws:states:::batch:submitJob.sync",
                "Parameters": {
                  "JobDefinition": "${var.job_definition_arn}",
                  "JobName": "cxg",
                  "JobQueue": "${var.job_queue_arn}",
                  "RetryStrategy": {
                    "Attempts": ${var.max_attempts},
                    "EvaluateOnExit": [
                      {
                        "Action": "EXIT",
                        "OnExitCode": "1"
                      },
                      {
                        "Action": "RETRY",
                        "OnExitCode": "*"
                      }
                    ]
                  },
                  "ContainerOverrides": {
                    "Environment": [
                      {
                        "Name": "DATASET_ID",
                        "Value.$": "$.dataset_id"
                      },
                      {
                        "Name": "STEP_NAME",
                        "Value": "cxg"
                      }
                    ]
                  }
                },
                "ResultPath": null,
                "TimeoutSeconds": 360000
              }
            }
          },
          {
            "StartAt": "Seurat",
            "States": {
              "Seurat": {
                "Type": "Task",
                "End": true,
                "Resource": "arn:aws:states:::batch:submitJob.sync",
                "Parameters": {
                  "JobDefinition": "${var.job_definition_arn}",
                  "JobName": "seurat",
                  "JobQueue": "${var.job_queue_arn}",
                  "RetryStrategy": {
                    "Attempts": ${var.max_attempts},
                    "EvaluateOnExit": [
                      {
                        "Action": "EXIT",
                        "OnExitCode": "1"
                      },
                      {
                        "Action": "RETRY",
                        "OnExitCode": "*"
                      }
                    ]
                  },
                  "ContainerOverrides": {
                    "Environment": [
                      {
                        "Name": "DATASET_ID",
                        "Value.$": "$.dataset_id"
                      },
                      {
                        "Name": "STEP_NAME",
                        "Value": "seurat"
                      }
                    ]
                  }
                },
                "TimeoutSeconds": 36000
              }
            }
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
      "HandleSuccess": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_success_handler}",
        "End": true,
        "Retry": [ {
            "ErrorEquals": ["Lambda.AWSLambdaException"],
            "IntervalSeconds": 1,
            "MaxAttempts": 3,
            "BackoffRate": 2.0
        } ]
      },
      "HandleErrors": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_error_handler}",
        "End": true,
        "Retry": [ {
            "ErrorEquals": ["Lambda.AWSLambdaException"],
            "IntervalSeconds": 1,
            "MaxAttempts": 3,
            "BackoffRate": 2.0
        } ]
      }
    }
}
EOF
}

resource "aws_sfn_state_machine" "state_machine_seurat" {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-seurat-sfn"
  role_arn = var.role_arn

  definition = <<EOF
{
  "StartAt": "Seurat",
  "States": {
    "Seurat": {
      "Type": "Task",
      "End": true,
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${var.job_definition_arn}",
        "JobName": "seurat",
        "JobQueue": "${var.job_queue_arn}",
        "ContainerOverrides": {
          "Environment": [
            {
              "Name": "DATASET_ID",
              "Value.$": "$.dataset_id"
            },
            {
              "Name": "STEP_NAME",
              "Value": "seurat"
            }
          ]
        }
      },
      "TimeoutSeconds": 36000
    }
  }
}
EOF
}

resource "aws_sfn_state_machine" "state_machine_cxg_remaster" {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-cxg-remaster-v2-sfn"
  role_arn = var.role_arn

  definition = <<EOF
{
  "StartAt": "CxgRemasterV2",
  "States": {
    "CxgRemasterV2": {
      "Type": "Task",
      "End": true,
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${var.job_definition_arn}",
        "JobName": "cxg_remaster",
        "JobQueue": "arn:aws:batch:us-west-2:699936264352:job-queue/dp-dev",
        "ContainerOverrides": {
          "Environment": [
            {
              "Name": "DATASET_ID",
              "Value.$": "$.dataset_id"
            },
            {
              "Name": "STEP_NAME",
              "Value": "cxg_remaster"
            }
          ]
        }
      },
      "TimeoutSeconds": 36000
    }
  }
}
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/upload-sfn"
}
