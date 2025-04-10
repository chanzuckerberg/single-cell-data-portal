output batch_job_definition {
  value       = aws_batch_job_definition.batch_job_def.id
  description = "ARN for the batch job definition"
}

output batch_job_definition_no_revision {
  value       = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/${aws_batch_job_definition.batch_job_def.name}"
  description = "ARN for the batch job definition"
}

output cxg_job_definition_no_revision {
  value       = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/${aws_batch_job_definition.cxg_job_def.name}"
  description = "ARN for the cxg batch job definition"
}

output atac_job_definition {
  value       = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/${aws_batch_job_definition.atac_job_def.name}"
  description = "ARN for the atac batch job definition"
}

output batch_job_log_group {
  value       = aws_cloudwatch_log_group.cloud_watch_logs_group.id
  description = "Name of the CloudWatch log group for the batch job"
}
