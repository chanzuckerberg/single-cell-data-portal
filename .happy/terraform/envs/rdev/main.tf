module stack {
  source                       = "./modules/ecs-stack"
  aws_account_id               = var.aws_account_id
  aws_role                     = var.aws_role
  happymeta_                   = var.happymeta_
  happy_config_secret          = var.happy_config_secret
  image_tag                    = var.image_tag
  priority                     = var.priority
  stack_name                   = var.stack_name
  deployment_stage             = "rdev"
  delete_protected             = false
  require_okta                 = true
  stack_prefix                 = "/${var.stack_name}"
  batch_container_memory_limit = 28000
  wmg_batch_container_memory_limit = 248000
  wmg_desired_vcpus                = 128
  cg_desired_vcpus                 = 128
  cg_batch_container_memory_limit  = 248000
  backend_instance_count       = 1
  frontend_instance_count      = 1
  backend_memory               = 6 * 1024
  backend_de_instance_count    = 1
  backend_de_memory            = 8192
  backend_de_cpus              = 1
  backend_de_workers           = 1
  backend_wmg_instance_count    = 1
  backend_wmg_memory            = 16 * 1024
  backend_wmg_cpus              = 4
  backend_wmg_workers           = 5  
  frontend_memory              = 4096
  dd_key_secret_arn            = "arn:aws:secretsmanager:us-west-2:699936264352:secret:dd_api_key-nGPNwx"

  wait_for_steady_state        = var.wait_for_steady_state
}
