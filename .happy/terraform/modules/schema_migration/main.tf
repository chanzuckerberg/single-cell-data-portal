data aws_region current {}

data aws_caller_identity current {}

locals {
  name = "schema-migration"
  job_definition_arn = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/dp-${var.deployment_stage}-${var.custom_stack_name}-schema-migration"
}

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
      {
        name= "DATASETS_BUCKET",
        value= var.datasets_bucket
      },
    ],
    resourceRequirements = [
        {
          type= "VCPU",
          Value="1"
        },
        {
          Type="MEMORY",
          Value = "8000"
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

resource aws_batch_job_definition dataset_migrations {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-dataset_migrations-${local.name}"
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
      {
        name= "DATASETS_BUCKET",
        value= var.datasets_bucket
      },
    ],
    resourceRequirements = [
        {
          type= "VCPU",
          Value="2"
        },
        {
          Type="MEMORY",
          Value = "16000"
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

resource aws_batch_job_definition pubish_revisions {
  type = "container"
  name = "dp-${var.deployment_stage}-${var.custom_stack_name}-${local.name}-publish-revisions"
  container_properties = jsonencode({
    command = ["python3",
      "-m",
      "backend.layers.processing.publish_revisions",
    ],
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
      {
        name= "DATASETS_BUCKET",
        value= var.datasets_bucket
      },
    ],
    resourceRequirements = [
      {
        type= "VCPU",
        Value="2"
      },
      {
        Type="MEMORY",
        Value = "4096"
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
          "limit_migration": "0"
        }
    },
    "ApplyDefaults": {
        "Type": "Pass",
        "Next": "GatherCollections",
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
              "Name": "LIMIT_MIGRATION",
              "Value.$": "$.limit_migration"
            }
          ]
        }
      },
      "Next": "SpanCollections",
      "Retry": [ {
          "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
          "IntervalSeconds": 2,
          "MaxAttempts": 7,
          "BackoffRate": 5
      }, {
          "ErrorEquals": ["States.TaskFailed"],
          "IntervalSeconds": 30,
          "MaxAttempts": 1
      }, {
          "ErrorEquals": ["States.ALL"],
          "IntervalSeconds": 2,
          "MaxAttempts": 3,
          "BackoffRate": 5
      }],
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
          "Mode": "DISTRIBUTED",
          "ExecutionType": "STANDARD"
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
                    "Name": "TASK_TOKEN",
                    "Value.$": "$$.Task.Token"
                  },
                  {
                    "Name": "EXECUTION_ID",
                    "Value.$": "$.execution_id"
                  }
                ]
              }
            },
            "Next": "DatasetsExists",
            "Retry": [ {
                "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
                "IntervalSeconds": 2,
                "MaxAttempts": 7,
                "BackoffRate": 5
            }, {
                "ErrorEquals": ["States.TaskFailed"],
                "IntervalSeconds": 30,
                "MaxAttempts": 1
            }, {
                "ErrorEquals": ["States.ALL"],
                "IntervalSeconds": 2,
                "MaxAttempts": 3,
                "BackoffRate": 5
            }],
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "ResultPath": null,
                "Next": "CollectionCleanup"
              }
            ]
          },
          "DatasetsExists": {
            "Type": "Choice",
            "Choices": [
              {
                "Variable": "$.key_name",
                "IsPresent": false,
                "Next": "CollectionCleanup"
              }
            ],
            "Default": "SpanDatasets"
          },
          "CollectionCleanup": {
            "Type": "Task",
            "Resource": "arn:aws:states:::batch:submitJob.sync",
            "Parameters": {
            "JobDefinition": "${resource.aws_batch_job_definition.schema_migrations.arn}",
            "JobName": "Collection_cleanup",
            "JobQueue": "${var.job_queue_arn}",
              "Timeout": {
                "AttemptDurationSeconds": 600
              },
              "ContainerOverrides": {
                "Environment": [
                  {
                    "Name": "STEP_NAME",
                    "Value": "collection_cleanup"
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
                    "Name": "TASK_TOKEN",
                    "Value.$": "$$.Task.Token"
                  },
                  {
                    "Name": "EXECUTION_ID",
                    "Value.$": "$.execution_id"
                  }
                ]
              }
            },
            "End": true,
            "Retry": [ {
                "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
                "IntervalSeconds": 2,
                "MaxAttempts": 7,
                "BackoffRate": 5
            }, {
                "ErrorEquals": ["States.TaskFailed"],
                "IntervalSeconds": 30,
                "MaxAttempts": 1
            }, {
                "ErrorEquals": ["States.ALL"],
                "IntervalSeconds": 2,
                "MaxAttempts": 3,
                "BackoffRate": 5
            }],
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "Next": "CollectionError",
                "ResultPath": null
              }
            ]
          },
          "SpanDatasets": {
            "Type": "Map",
            "ItemProcessor": {
              "ProcessorConfig": {
                "Mode": "DISTRIBUTED",
                "ExecutionType": "STANDARD"
              },
              "StartAt": "DatasetMigration",
              "States": {
                "DatasetMigration": {
                  "Type": "Task",
                  "Resource": "arn:aws:states:::batch:submitJob.sync",
                  "Parameters": {
                    "JobDefinition": "${resource.aws_batch_job_definition.dataset_migrations.arn}",
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
                          "Value.$": "$.execution_id"
                        }
                      ]
                    }
                  },
                  "Next": "StepFunctionsStartExecution",
                  "Retry": [ {
                      "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
                      "IntervalSeconds": 2,
                      "MaxAttempts": 7,
                      "BackoffRate": 5
                  }, {
                      "ErrorEquals": ["States.TaskFailed"],
                      "IntervalSeconds": 30,
                      "MaxAttempts": 1
                  }, {
                      "ErrorEquals": ["States.ALL"],
                      "IntervalSeconds": 2,
                      "MaxAttempts": 3,
                      "BackoffRate": 5
                  }],
                  "Catch": [
                    {
                      "ErrorEquals": [
                        "States.ALL"
                      ],
                      "Next": "DatasetError",
                      "ResultPath": null
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
                      "manifest.$": "$.result.manifest",
                      "dataset_version_id.$": "$.result.dataset_version_id",
                      "collection_version_id.$": "$.result.collection_version_id",
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
                      "ResultPath": null
                    }
                  ],
                  "ResultPath": null
                }
              }
            },
            "ItemReader": {
              "Resource": "arn:aws:states:::s3:getObject",
              "ReaderConfig": {
                "InputType": "JSON"
              },
              "Parameters": {
                "Bucket": "${var.artifact_bucket}",
                "Key.$": "$.key_name"
              }
            },
            "Next": "CollectionCleanup",
            "MaxConcurrency": 10,
            "Catch": [
              {
                "ErrorEquals": [
                  "States.ALL"
                ],
                "Next": "CollectionCleanup",
                "ResultPath": null
              }
            ],
            "OutputPath": "$[0]",
            "ToleratedFailurePercentage": 100
          },
          "CollectionError": {
            "Type": "Pass",
            "End": true
          }
        }
      },
      "ItemReader": {
        "Resource": "arn:aws:states:::s3:getObject",
        "ReaderConfig": {
          "InputType": "JSON"
        },
        "Parameters": {
          "Bucket": "${var.artifact_bucket}",
          "Key.$": "$.key_name"
        }
      },
      "MaxConcurrency": 10,
      "Next": "report",
      "Catch": [
        {
          "ErrorEquals": [
            "States.ALL"
          ],
          "Next": "report",
          "ResultPath": null
        }
      ],
      "ToleratedFailurePercentage": 20
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
      "Retry": [ {
          "ErrorEquals": ["AWS.Batch.TooManyRequestsException", "Batch.BatchException", "Batch.AWSBatchException"],
          "IntervalSeconds": 2,
          "MaxAttempts": 7,
          "BackoffRate": 5
      }, {
          "ErrorEquals": ["States.TaskFailed"],
          "IntervalSeconds": 30,
          "MaxAttempts": 1
      }, {
          "ErrorEquals": ["States.ALL"],
          "IntervalSeconds": 2,
          "MaxAttempts": 3,
          "BackoffRate": 5
      }],
      "End": true
    }
  }
}
EOF
}
