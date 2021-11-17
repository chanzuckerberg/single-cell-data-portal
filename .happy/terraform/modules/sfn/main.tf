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
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "DROPBOX_URL",
                "Value.$": "$.url"
              },
              {
                "Name": "DATASET_ID",
                "Value.$": "$.dataset_uuid"
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
      "CxgSeuratParallel": {
        "Type": "Parallel",
        "End": true,
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
                  "ContainerOverrides": {
                    "Environment": [
                      {
                        "Name": "DATASET_ID",
                        "Value.$": "$.dataset_uuid"
                      },
                      {
                        "Name": "STEP_NAME",
                        "Value": "cxg"
                      }
                    ]
                  }
                },
                "ResultPath": null,
                "TimeoutSeconds": 36000,
                "Retry": [
                  {
                    "ErrorEquals": [
                      "States.TaskFailed"
                    ],
                    "IntervalSeconds": 1,
                    "BackoffRate": 2,
                    "MaxAttempts": 2
                  }
                ]
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
                  "JobName": "cxg",
                  "JobQueue": "${var.job_queue_arn}",
                  "ContainerOverrides": {
                    "Environment": [
                      {
                        "Name": "DATASET_ID",
                        "Value.$": "$.dataset_uuid"
                      },
                      {
                        "Name": "STEP_NAME",
                        "Value": "seurat"
                      }
                    ]
                  }
                },
                "TimeoutSeconds": 36000,
                "Retry": [
                  {
                    "ErrorEquals": [
                      "States.TaskFailed"
                    ],
                    "IntervalSeconds": 1,
                    "BackoffRate": 2,
                    "MaxAttempts": 2
                  }
                ]
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
        "JobName": "cxg",
        "JobQueue": "${var.job_queue_arn}",
        "ContainerOverrides": {
          "Environment": [
            {
              "Name": "DATASET_ID",
              "Value.$": "$.dataset_uuid"
            },
            {
              "Name": "STEP_NAME",
              "Value": "seurat"
            }
          ]
        }
      },
      "TimeoutSeconds": 36000,
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
