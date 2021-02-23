output batch_job_definition {
  value       = aws_batch_job_definition.batch_job_def.id
  description = "ARN for the batch job definition"
}
