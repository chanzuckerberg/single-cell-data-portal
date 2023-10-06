data aws_region current {}

data aws_caller_identity current {}

locals {
  name = "schema-migration"
  job_definition_arn = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/dp-${var.deployment_stage}-${var.custom_stack_name}-schema-migration"
  swap_job_definition_arn = "${local.job_definition_arn}-swap"
}

resource aws_cloudwatch_log_group batch_cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/${local.name}-batch"
}

resource aws_batch_job_definition schema_migrations_swap {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-${local.name}-swap"
  container_properties = jsonencode({
    jobRoleArn= var.batch_role_arn,
    image= var.image,
    environment= [
      {
        name= "ARTIFACT_BUCKET",
        value= var.artifact_bucket
      },
      {
        name= "DEPLOYMENT_STAGE",
        value= var.deployment_stage
      },
      {
        name= "AWS_DEFAULT_REGION",
        value= data.aws_region.current.name
      },
      {
        name= "REMOTE_DEV_PREFIX",
        value= var.remote_dev_prefix
      },
    ],
    resourceRequirements = [
      {
        type= "VCPU",
        Value="6"
      },
      {
        Type="MEMORY",
        Value = "47500"
      }
    ]
    linuxParameters= {
     maxSwap= 200000,
     swappiness= 60
    },
    logConfiguration= {
      logDriver= "awslogs",
      options= {
        awslogs-group= aws_cloudwatch_log_group.batch_cloud_watch_logs_group.id,
        awslogs-region= data.aws_region.current.name
      }
    }
  })
}

resource aws_batch_job_definition schema_migrations {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-${local.name}"
  container_properties = jsonencode({
    jobRoleArn= var.batch_role_arn,
    image= var.image,
    environment= [
      {
        name= "ARTIFACT_BUCKET",
        value= var.artifact_bucket
      },
      {
        name= "DEPLOYMENT_STAGE",
        value= var.deployment_stage
      },
      {
        name= "AWS_DEFAULT_REGION",
        value= data.aws_region.current.name
      },
      {
        name= "REMOTE_DEV_PREFIX",
        value= var.remote_dev_prefix
      },
    ],
    resourceRequirements = [
                  {
              type= "VCPU",
              Value="2"
            },
            {
              Type="MEMORY",
              Value = "2048"
            }
    ]
    logConfiguration= {
      logDriver= "awslogs",
      options= {
        awslogs-group= aws_cloudwatch_log_group.batch_cloud_watch_logs_group.id,
        awslogs-region= data.aws_region.current.name
      }
    }
  })
}

resource aws_sfn_state_machine sfn_schema_migration {
  name     = "dp-${var.deployment_stage}-${var.custom_stack_name}-${local.name}-sfn"
  role_arn = var.sfn_role_arn

  definition = <<EOF
{
  "Comment": "Schema Migration State Machine",
  "StartAt": "DefineDefaults",
  "States": {
    "DefineDefaults": {
        "Type": "Pass",
        "Next": "ApplyDefaults",
        "ResultPath": "$.inputDefaults",
        "Parameters": {
          "auto_publish": "False"
        }
    },
    "ApplyDefaults": {
        "Type": "Pass",
        "Next": "DownloadValidate",
        "Parameters": {
          "args.$": "States.JsonMerge($.inputDefaults, $$.Execution.Input, false)"
        },
        "ResultPath": "$.withDefaults",
        "OutputPath": "$.withDefaults.args"
    },
    "GatherCollections": {
      "Type": "Task",
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}",
        "JobName": "gather_collections",
        "JobQueue": "${var.job_queue_arn}",
        "Timeout": {
          "AttemptDurationSeconds": 600
        },
        "ContainerOverrides": {
          "Environment": [
            {
              "Name": "STEP_NAME",
              "Value": "gather_collections"
            },
            {
              "Name": "MIGRATE",
              "Value": "True"
            },
            {
              "Name": "TASK_TOKEN",
              "Value.$": "$$.Task.Token"
            },
            {
              "Name": "EXECUTION_ID",
              "Value.$": "$$.Execution.Name"
            },
            {
              "Name": "AUTO_PUBLISH",
              "Value.$": "$.auto_publish"
            }
          ]
        }
      },
      "Next": "SpanCollections",
      "Catch": [
        {
          "ErrorEquals": [
            "States.ALL"
          ],
          "Next": "report",
          "ResultPath": null
        }
      ]
    },
    "SpanCollections": {
      "Type": "Map",
      "ItemProcessor": {
        "ProcessorConfig": {
          "Mode": "INLINE"
        },
        "StartAt": "CollectionMigration",
        "States": {
          "CollectionMigration": {
            "Type": "Task",
            "Resource": "arn:aws:states:::batch:submitJob.sync",
            "Parameters": {
              "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}",
              "JobName": "collection_migration",
              "JobQueue": "${var.job_queue_arn}",
              "Timeout": {
                "AttemptDurationSeconds": 600
              },
              "ContainerOverrides": {
                "Environment": [
                  {
                    "Name": "STEP_NAME",
                    "Value": "collection_migrate"
                  },
                  {
                    "Name": "MIGRATE",
                    "Value": "True"
                  },
                  {
                    "Name": "COLLECTION_ID",
                    "Value.$": "$.collection_id"
                  },
                  {
                    "Name": "COLLECTION_VERSION_ID",
                    "Value.$": "$.collection_version_id"
                  },
                  {
                    "Name": "CAN_PUBLISH",
                    "Value.$": "$.can_publish"
                  },
                  {
                    "Name": "TASK_TOKEN",
                    "Value.$": "$$.Task.Token"
                  },
                  {
                    "Name": "EXECUTION_ID",
                    "Value.$": "$$.Execution.Name"
                  }
                ]
              }
            },
            "Next": "DatasetsExists",
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "ResultPath": "$.error",
                "Next": "CollectionPublish"
              }
            ]
          },
          "DatasetsExists": {
            "Type": "Choice",
            "Choices": [
              {
                "Variable": "$.no_datasets",
                "IsPresent": true,
                "Next": "CollectionPublish"
              }
            ],
            "Default": "SpanDatasets"
          },
          "CollectionPublish": {
            "Type": "Task",
            "Resource": "arn:aws:states:::batch:submitJob.sync",
            "Parameters": {
            "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}",
            "JobName": "Collection_publish",
            "JobQueue": "${var.job_queue_arn}",
              "Timeout": {
                "AttemptDurationSeconds": 600
              },
              "ContainerOverrides": {
                "Environment": [
                  {
                    "Name": "STEP_NAME",
                    "Value": "collection_publish"
                  },
                  {
                    "Name": "MIGRATE",
                    "Value": "True"
                  },
                  {
                    "Name": "COLLECTION_VERSION_ID",
                    "Value.$": "$.collection_version_id"
                  },
                  {
                    "Name": "CAN_PUBLISH",
                    "Value.$": "$.can_publish"
                  },
                  {
                    "Name": "TASK_TOKEN",
                    "Value.$": "$$.Task.Token"
                  },
                  {
                    "Name": "EXECUTION_ID",
                    "Value.$": "$$.Execution.Name"
                  }
                ]
              }
            },
            "End": true,
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "Next": "CollectionError",
                "ResultPath": "$.error"
              }
            ]
          },
          "SpanDatasets": {
            "Type": "Map",
            "ItemProcessor": {
              "ProcessorConfig": {
                "Mode": "INLINE"
              },
              "StartAt": "DatasetMigration",
              "States": {
                "DatasetMigration": {
                  "Type": "Task",
                  "Resource": "arn:aws:states:::batch:submitJob.sync",
                  "Parameters": {
                    "JobDefinition": "${local.swap_job_definition_arn}",
                    "JobName": "dataset_migration",
                    "JobQueue": "${var.job_queue_arn}",
                    "Timeout": {
                      "AttemptDurationSeconds": 36000
                    },
                    "ContainerOverrides": {
                      "Environment": [
                        {
                          "Name": "STEP_NAME",
                          "Value": "dataset_migrate"
                        },
                        {
                          "Name": "MIGRATE",
                          "Value": "True"
                        },
                        {
                          "Name": "DATASET_ID",
                          "Value.$": "$.dataset_id"
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
                          "Name": "COLLECTION_ID",
                          "Value.$": "$.collection_id"
                        },
                        {
                          "Name": "TASK_TOKEN",
                          "Value.$": "$$.Task.Token"
                        },
                        {
                          "Name": "EXECUTION_ID",
                          "Value.$": "$$.Execution.Name"
                        }
                      ]
                    }
                  },
                  "Next": "StepFunctionsStartExecution",
                  "Catch": [
                    {
                      "ErrorEquals": [
                        "States.ALL"
                      ],
                      "Next": "DatasetError",
                      "ResultPath": "$.error"
                    }
                  ],
                  "ResultPath": "$.result"
                },
                "DatasetError": {
                  "Type": "Pass",
                  "End": true
                },
                "StepFunctionsStartExecution": {
                  "Type": "Task",
                  "Resource": "arn:aws:states:::states:startExecution.sync:2",
                  "Parameters": {
                    "StateMachineArn": "arn:aws:states:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:stateMachine:dp-${var.deployment_stage}-${var.custom_stack_name}-sfn",
                    "Name.$": "$.result.sfn_name",
                    "Input": {
                      "AWS_STEP_FUNCTIONS_STARTED_BY_EXECUTION_ID.$": "$$.Execution.Id",
                      "url.$": "$.result.uri",
                      "dataset_id.$": "$.result.dataset_version_id",
                      "collection_id.$": "$.result.collection_version_id",
                      "job_queue": "${var.job_queue_arn}"
                    }
                  },
                  "End": true,
                  "Catch": [
                    {
                      "ErrorEquals": [
                        "States.ALL"
                      ],
                      "Next": "DatasetError",
                      "ResultPath": "$.error"
                    }
                  ],
                  "ResultPath": "$.result"
                }
              }
            },
            "ItemsPath": "$.datasets",
            "Next": "CollectionPublish",
            "MaxConcurrency": 32,
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "Next": "CollectionPublish",
                "ResultPath": "$.error"
              }
            ],
            "OutputPath": "$[0]"
          },
          "CollectionError": {
            "Type": "Pass",
            "End": true
          }
        }
      },
      "ItemsPath": "$",
      "MaxConcurrency": 40,
      "Next": "report",
      "Catch": [
        {
          "ErrorEquals": [
            "States.ALL"
          ],
          "Next": "report",
          "ResultPath": null
        }
      ]
    },
    "report": {
      "Type": "Task",
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}",
        "JobName": "report",
        "JobQueue": "${var.job_queue_arn}",
        "Timeout": {
          "AttemptDurationSeconds": 600
        },
        "ContainerOverrides": {
          "ResourceRequirements": [
            {
              "Type": "MEMORY",
              "Value": "8048"
            }
          ],
          "Environment": [
            {
              "Name": "STEP_NAME",
              "Value": "report"
            },
            {
              "Name": "MIGRATE",
              "Value": "True"
            },
            {
              "Name": "TASK_TOKEN",
              "Value.$": "$$.Task.Token"
            },
            {
              "Name": "EXECUTION_ID",
              "Value.$": "$$.Execution.Name"
            }
          ]
        }
      },
      "End": true
    }
  }
}
EOF
}
