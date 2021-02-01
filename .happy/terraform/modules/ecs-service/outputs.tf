output frontend_url {
  value       = module.dev_env.frontend_url
  description = "The URL endpoint for the website service"
}

output backend_url {
  value       = module.dev_env.backend_url
  description = "The URL endpoint for the website service"
}

output delete_db_task_definition_arn {
  value       = module.dev_env.delete_db_task_definition_arn
  description = "ARN of the Deletion ECS Task Definition"
}

output migrate_db_task_definition_arn {
  value       = module.dev_env.migrate_db_task_definition_arn
  description = "ARN of the Migration ECS Task Definition"
}
