# This is used for environment (dev, staging, prod) deployments
locals {
  h5ad_timeout = 86400 # 24 hours
  cxg_timeout = 172800 # 48 hours
}

data aws_region current {}

resource "aws_sfn_state_machine" "state_machine" {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-sfn"
  role_arn = var.role_arn

  definition = <<EOF
{
    "StartAt": "DefineDefaults",
    "States": {
      "DefineDefaults": {
        "Type": "Pass",
        "Next": "ApplyDefaults",
        "ResultPath": "$.inputDefaults",
        "Parameters": {
          "job_queue": "${var.job_queue_arn}"
        }
      },
      "ApplyDefaults": {
        "Type": "Pass",
        "Next": "ValidateAnndata",
        "Parameters": {
          "args.$": "States.JsonMerge($.inputDefaults, $$.Execution.Input, false)"
        },
        "ResultPath": "$.withDefaults",
        "OutputPath": "$.withDefaults.args"
      },
      "ValidateAnndata": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Next": "AddLabels",
        "Parameters": {
          "JobDefinition":"${var.job_definition_arn}",
          "JobName": "validate_anndata",
          "JobQueue.$": "$.job_queue",
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "MANIFEST",
                "Value.$": "$.manifest"
              },
              {
                "Name": "DATASET_VERSION_ID",
                "Value.$": "$.dataset_version_id"
              },
              {
                "Name": "COLLECTION_VERSION_ID",
                "Value.$": "$.collection_version_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "validate_anndata"
              }
            ]
          }
        },
        "ResultPath": null,
        "TimeoutSeconds": ${local.h5ad_timeout},
        "Retry": [ {
            "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
            "IntervalSeconds": 2,
            "MaxAttempts": 7,
            "BackoffRate": 5
        } ],
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
      "AddLabels": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Next": "Cxg",
        "Parameters": {
          "JobDefinition": "${var.job_definition_arn}",
          "JobName": "add_labels",
          "JobQueue.$": "$.job_queue",
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "DATASET_VERSION_ID",
                "Value.$": "$.dataset_version_id"
              },
              {
                "Name": "COLLECTION_VERSION_ID",
                "Value.$": "$.collection_version_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "add_labels"
              }
            ]
          }
        },
        "ResultPath": null,
        "TimeoutSeconds": ${local.h5ad_timeout},
        "Retry": [ {
            "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
            "IntervalSeconds": 2,
            "MaxAttempts": 7,
            "BackoffRate": 5
        } ],
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
      "Cxg": {
        "Type": "Task",
        "Next": "HandleSuccess",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Parameters": {
          "JobDefinition":"${var.cxg_definition_arn}",
          "JobName": "cxg",
          "JobQueue.$": "$.job_queue",
          "ContainerOverrides": {
            "Environment": [
              {
                "Name": "DATASET_VERSION_ID",
                "Value.$": "$.dataset_version_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "cxg"
              }
            ]
          }
        },
        "Retry": [ {
            "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
            "IntervalSeconds": 2,
            "MaxAttempts": 7,
            "BackoffRate": 5
        } ],
        "Catch": [
          {
            "ErrorEquals": [
              "States.ALL"
            ],
            "Next": "HandleErrors",
            "ResultPath": "$.error"
          }
        ],
        "ResultPath": null,
        "TimeoutSeconds": ${local.cxg_timeout}
      },
      "HandleSuccess": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_success_handler}",
        "Parameters": {
          "execution_id.$": "$$.Execution.Id",
          "cxg_job.$": "$"
        },
        "Retry": [ {
            "ErrorEquals": ["Lambda.AWSLambdaException"],
            "IntervalSeconds": 1,
            "MaxAttempts": 3,
            "BackoffRate": 2.0
        } ],
        "End": true,
        "ResultPath": null
      },
      "HandleErrors": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_error_handler}",
        "Parameters": {
          "execution_id.$": "$$.Execution.Id",
          "error.$": "$.error",
          "dataset_version_id.$": "$.dataset_version_id",
          "collection_version_id.$": "$.collection_version_id"
        },
        "Retry": [ {
            "ErrorEquals": ["Lambda.AWSLambdaException"],
            "IntervalSeconds": 1,
            "MaxAttempts": 3,
            "BackoffRate": 2.0
        } ],
        "ResultPath": null,
        "Next": "RaiseError"
      },
      "RaiseError": {
        "Type": "Fail",
        "Cause": "Failed to ingest dataset.",
        "End": true
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
        "JobDefinition": "${var.cxg_definition_arn}",
        "JobName": "cxg_remaster",
        "JobQueue": "${var.job_queue_arn}",
        "ContainerOverrides": {
          "Environment": [
            {
              "Name": "DATASET_VERSION_ID",
              "Value.$": "$.dataset_version_id"
            },
            {
              "Name": "STEP_NAME",
              "Value": "cxg_remaster"
            }
          ]
        }
      },
      "TimeoutSeconds": ${local.cxg_timeout}
    }
  }
}
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/upload-sfn"
}
