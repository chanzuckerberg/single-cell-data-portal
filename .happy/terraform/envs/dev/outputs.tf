output frontend_url {
  value       = module.stack.frontend_url
  description = "The URL endpoint for the website service"
}

output backend_url {
  value       = module.stack.backend_url
  description = "The URL endpoint for the website service"
}

output delete_db_task_definition_arn {
  value       = module.stack.delete_db_task_definition_arn
  description = "ARN of the Deletion ECS Task Definition"
}

output migrate_db_task_definition_arn {
  value       = module.stack.migrate_db_task_definition_arn
  description = "ARN of the Migration ECS Task Definition"
}
