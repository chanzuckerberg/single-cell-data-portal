locals {
  name = "schema-migration"
}

data aws_region current {}

data aws_caller_identity current {}

resource aws_cloudwatch_log_group batch_cloud_watch_logs_group {
  retention_in_days = 365
  name              = "/dp/${var.deployment_stage}/${var.custom_stack_name}/${local.name}-batch"
}

resource aws_batch_job_definition schema_migrations {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-${local.name}"
  container_properties = jsonencode({
    jobRoleArn= var.batch_role_arn,
    image= var.image,
    memory= var.batch_container_memory_limit,
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
    vcpus= 2,
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
  "StartAt": "gather collections",
  "States": {
    "gather collections": {
      "Type": "Task",
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}:$LATEST",
        "JobName": "gather_collections",
        "JobQueue": "${var.job_queue_arn}",
        "Timeout": {
          "AttemptDurationSeconds": 600
        },
        "ContainerOverrides": {
          "ResourceRequirements": [
            {
              "Type": "VCPU",
              "Value": "4"
            },
            {
              "Type": "MEMORY",
              "Value": "2048"
            }
          ],
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
            }
          ]
        }
      },
      "Next": "span collections"
    },
    "span collections": {
      "Type": "Map",
      "ItemProcessor": {
        "ProcessorConfig": {
          "Mode": "INLINE"
        },
        "StartAt": "collection migration",
        "States": {
          "collection migration": {
            "Type": "Task",
            "Resource": "arn:aws:states:::batch:submitJob.sync",
            "Parameters": {
              "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}:$LATEST",
              "JobName": "collection_migration",
              "JobQueue": "${var.job_queue_arn}",
              "Timeout": {
                "AttemptDurationSeconds": 600
              },
              "ContainerOverrides": {
                "ResourceRequirements": [
                  {
                    "Type": "VCPU",
                    "Value": "4"
                  },
                  {
                    "Type": "MEMORY",
                    "Value": "2048"
                  }
                ],
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
                    "Name": "CAN_OPEN_REVISION",
                    "Value.$": "$.can_open_revision"
                  },
                  {
                      "Name": "TASK_TOKEN",
                      "Value.$": "$$.Task.Token"
                  }
                ]
              }
            },
            "Next": "span datasets"
          },
          "span datasets": {
            "Type": "Map",
            "ItemProcessor": {
              "ProcessorConfig": {
                "Mode": "INLINE"
              },
              "StartAt": "dataset_migration",
              "States": {
                "dataset_migration": {
                  "Type": "Task",
                  "Resource": "arn:aws:states:::batch:submitJob.sync",
                  "Parameters": {
                    "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}:$LATEST",
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
                          "Name": "COLLECTION_VERSION_ID",
                          "Value.$": "$.collection_version_id"
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
                            "Name": "TASK_TOKEN",
                            "Value.$": "$$.Task.Token"
                        }
                      ]
                    }
                  },
                  "Next": "Step Functions StartExecution"
                },
                "Step Functions StartExecution": {
                  "Type": "Task",
                  "Resource": "arn:aws:states:::states:startExecution.sync:2",
                  "Parameters": {
                    "StateMachineArn": "arn:aws:states:REGION:ACCOUNT_ID:stateMachine:STATE_MACHINE_NAME",
                    "Input": {
                      "AWS_STEP_FUNCTIONS_STARTED_BY_EXECUTION_ID.$": "$$.Execution.Id",
                      "url": "$.url",
                      "dataset_id": "$.dataset_version_id",
                      "collection_id": "$.collection_id",
                    }
                  },
                  "End": true
                }
              }
            },
            "ItemsPath": "$.datasets",
            "Next": "collection publish",
            "MaxConcurrency": 40
          },
          "collection publish": {
            "Type": "Task",
            "Resource": "arn:aws:states:::batch:submitJob",
            "Parameters": {
              "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}:$LATEST",
              "JobName": "collection_publish",
              "JobQueue": "${var.job_queue_arn}",
              "Timeout": {
                "AttemptDurationSeconds": 600
              },
              "ContainerOverrides": {
                "ResourceRequirements": [
                  {
                    "Type": "VCPU",
                    "Value": "4"
                  },
                  {
                    "Type": "MEMORY",
                    "Value": "2048"
                  }
                ],
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
                    "Name": "COLLECTION_ID",
                    "Value": "$.collection_id"
                  },
                  {
                    "Name": "CAN_PUBLISH",
                    "Value": "$.can_publish"
                  },
                  {
                      "Name": "TASK_TOKEN",
                      "Value.$": "$$.Task.Token"
                  }
                ]
              }
            },
            "End": true
          }
        }
      },
      "ItemsPath": "$.collections",
      "MaxConcurrency": 40,
      "Next": "report"
    },
    "report": {
      "Type": "Task",
      "Resource": "arn:aws:states:::batch:submitJob.sync",
      "Parameters": {
        "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}:$LATEST",
        "JobName": "report",
        "JobQueue": "${var.job_queue_arn}",
        "Timeout": {
          "AttemptDurationSeconds": 600
        },
        "ContainerOverrides": {
          "ResourceRequirements": [
            {
              "Type": "VCPU",
              "Value": "4"
            },
            {
              "Type": "MEMORY",
              "Value": "2048"
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
