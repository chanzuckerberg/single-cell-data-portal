output task_definition_arn {
  value       = aws_ecs_task_definition.task_definition.arn
  description = "ARN of the Migration ECS Task Definition"
}
