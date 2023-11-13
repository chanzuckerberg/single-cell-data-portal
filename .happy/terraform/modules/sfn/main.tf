# Same file as https://github.com/chanzuckerberg/single-cell-infra/blob/main/.happy/terraform/modules/sfn/main.tf
# This is used for environment (dev, staging, prod) deployments
locals {
  timeout = 86400 # 24 hours
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
        "Next": "Download",
        "Parameters": {
          "args.$": "States.JsonMerge($.inputDefaults, $$.Execution.Input, false)"
        },
        "ResultPath": "$.withDefaults",
        "OutputPath": "$.withDefaults.args"
      },
      "Download": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Next": "RegisterJobDefinition",
        "Parameters": {
          "JobDefinition":"${var.job_definition_arn}",
          "JobName": "download",
          "JobQueue.$": "$.job_queue",
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
                "Name": "DATASET_VERSION_ID",
                "Value.$": "$.dataset_version_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "download"
              },
              {
                "Name": "TASK_TOKEN",
                "Value.$": "$$.Task.Token"
              }
            ]
          }
        },
        "TimeoutSeconds": ${local.timeout},
        "Catch": [
          {
            "ErrorEquals": [
              "States.ALL"
            ],
            "Next": "HandleErrors",
            "ResultPath": "$.error"
          }
        ],
        "ResultPath": "$.batch"
      },
      "RegisterJobDefinition": {
        "Type": "Task",
        "Next": "Validate",
        "Parameters": {
          "JobDefinitionName.$": "$.batch.JobDefinitionName",
          "Type": "container",
          "ContainerProperties" :{
            "Image" : "${var.image}",
            "JobRoleArn": "${var.batch_role_arn}",
            "Environment" : [
              {
                "Name" : "ARTIFACT_BUCKET",
                "Value" : "${var.artifact_bucket}"
              },
              {
                "Name" : "CELLXGENE_BUCKET",
                "Value" : "${var.cellxgene_bucket}"
              },
              {
                "Name" : "DATASETS_BUCKET",
                "Value" : "${var.datasets_bucket}"
              },
              {
                "Name" : "DEPLOYMENT_STAGE",
                "Value" : "${var.deployment_stage}"
              },
              {
                "Name" : "AWS_DEFAULT_REGION",
                "Value" : "${data.aws_region.current.name}"
              },
              {
                "Name" : "REMOTE_DEV_PREFIX",
                "Value" : "${var.remote_dev_prefix}"
              },
              {
                "Name" : "FRONTEND_URL",
                "Value" : "${var.frontend_url}"
              }
            ],
            "Vcpus.$" : "$.batch.Vcpus",
            "Memory.$" : "$.batch.Memory",
            "LinuxParameters.$" : "$.batch.LinuxParameters",
            "LogConfiguration" : {
              "LogDriver" : "awslogs",
              "Options" : {
                "awslogs-group" : "${var.batch_job_log_group}",
                "awslogs-region" : "${data.aws_region.current.name}"
              }
            }
          }
        },
        "Resource": "arn:aws:states:::aws-sdk:batch:registerJobDefinition",
        "ResultPath": "$.batch"
      },
      "Validate": {
        "Type": "Task",
        "Resource": "arn:aws:states:::batch:submitJob.sync",
        "Next": "CxgSeuratParallel",
        "Parameters": {
          "JobDefinition.$": "$.batch.JobDefinitionName",
          "JobName": "validate",
          "JobQueue.$": "$.job_queue",
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
                "Name": "DATASET_VERSION_ID",
                "Value.$": "$.dataset_version_id"
              },
              {
                "Name": "COLLECTION_VERSION_ID",
                "Value.$": "$.collection_version_id"
              },
              {
                "Name": "STEP_NAME",
                "Value": "validate"
              }
            ]
          }
        },
        "ResultPath": null,
        "TimeoutSeconds": ${local.timeout},
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
        "Branches": [
          {
            "StartAt": "Cxg",
            "States": {
              "Cxg": {
                "Type": "Task",
                "End": true,
                "Resource": "arn:aws:states:::batch:submitJob.sync",
                "Parameters": {
                  "JobDefinition.$": "$.batch.JobDefinitionName",
                  "JobName": "cxg",
                  "JobQueue.$": "$.job_queue",
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
                "Catch": [
                  {
                    "ErrorEquals": [
                      "States.ALL"
                    ],
                    "Next": "CatchCxgFailure",
                    "ResultPath": "$.error"
                  }
                ],
                "ResultPath": null,
                "TimeoutSeconds": 360000
              },
              "CatchCxgFailure": {
                "Type": "Pass",
                "End": true
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
                  "JobDefinition.$": "$.batch.JobDefinitionName",
                  "JobName": "seurat",
                  "JobQueue.$": "$.job_queue",
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
                        "Name": "DATASET_VERSION_ID",
                        "Value.$": "$.dataset_version_id"
                      },
                      {
                        "Name": "STEP_NAME",
                        "Value": "seurat"
                      }
                    ]
                  }
                },
                "Catch": [
                  {
                    "ErrorEquals": [
                      "States.ALL"
                    ],
                    "Next": "CatchSeuratFailure",
                    "ResultPath": "$.error"
                  }
                ],
                "TimeoutSeconds": ${local.timeout}
              },
              "CatchSeuratFailure": {
                "Type": "Pass",
                "End": true
              }
            }
          }
        ]
      },
      "HandleSuccess": {
        "Type": "Task",
        "InputPath": "$",
        "Resource": "${var.lambda_success_handler}",
        "Parameters": {
          "execution_id.$": "$$.Execution.Id",
          "cxg_job.$": "$[0]",
          "seurat_job.$": "$[1]"
        },
        "Retry": [ {
            "ErrorEquals": ["Lambda.AWSLambdaException"],
            "IntervalSeconds": 1,
            "MaxAttempts": 3,
            "BackoffRate": 2.0
        } ],
        "Next": "DeregisterJobDefinition",
        "ResultPath": null,
        "OutputPath": "$.[0]"
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
        "Next": "DeregisterJobDefinition",
        "ResultPath": null
      },
      "DeregisterJobDefinition": {
        "Type": "Task",
        "End": true,
        "Parameters": {
          "JobDefinition.$": "$.batch.JobDefinitionName"
        },
        "Resource": "arn:aws:states:::aws-sdk:batch:deregisterJobDefinition"
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
              "Name": "DATASET_VERSION_ID",
              "Value.$": "$.dataset_version_id"
            },
            {
              "Name": "STEP_NAME",
              "Value": "seurat"
            }
          ]
        }
      },
      "TimeoutSeconds": ${local.timeout}
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
      "TimeoutSeconds": ${local.timeout}
    }
  }
}
EOF
}

resource aws_cloudwatch_log_group cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/upload-sfn"
}
