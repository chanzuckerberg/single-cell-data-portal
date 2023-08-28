output batch_job_definition {
  value       = aws_batch_job_definition.batch_job_def.id
  description = "ARN for the CellGuide batch job definition"
}

output batch_job_definition_no_revision {
  value       = "arn:aws:batch:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:job-definition/${aws_batch_job_definition.batch_job_def.name}"
  description = "ARN for the CellGuide batch job definition"
}