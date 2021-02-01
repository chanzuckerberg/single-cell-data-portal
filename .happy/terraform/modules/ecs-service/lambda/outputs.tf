output error_handler {
  value       = aws_lambda_function.lambda_job_def.arn
  description = "ARN for the lambda error handler"
}
