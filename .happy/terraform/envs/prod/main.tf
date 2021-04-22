module stack {
  source                       = "./modules/ecs-stack"
  aws_account_id               = var.aws_account_id
  aws_role                     = var.aws_role
  happymeta_                   = var.happymeta_
  happy_config_secret          = var.happy_config_secret
  image_tag                    = var.image_tag
  priority                     = var.priority
  stack_name                   = var.stack_name
  deployment_stage             = "prod"
  delete_protected             = true
  require_okta                 = false
  frontend_url                 = "https://cellxgene.cziscience.com"
  backend_url                  = "https://api.cellxgene.cziscience.com"
  stack_prefix                 = ""
  batch_container_memory_limit = 250000
  backend_instance_count       = 6

  wait_for_steady_state        = var.wait_for_steady_state
}
