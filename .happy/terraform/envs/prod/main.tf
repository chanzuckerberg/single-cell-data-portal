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
  batch_container_memory_limit = 63500
  wmg_batch_container_memory_limit = 248000
  wmg_desired_vcpus                = 128
  cg_desired_vcpus                 = 128
  cg_batch_container_memory_limit  = 248000
  backend_memory               = 6 * 1024
  frontend_memory              = 4 * 1024
  backend_instance_count       = 2
  backend_cpus                 = 2
  backend_workers              = 8 # Rule of thumb we are using is 2*num CPUs+3 since backend is I/O bound
  backend_de_instance_count    = 3
  backend_de_memory            = 16 * 1024
  backend_de_cpus              = 4 # Rule of thumb we are using is num CPUs+1 since backend_de is compute bound
  backend_de_workers           = 5
  backend_wmg_instance_count    = 3
  backend_wmg_memory            = 16 * 1024
  backend_wmg_cpus              = 3
  backend_wmg_workers           = 7 # Rule of thumb we are using is num 2*CPUs+1 since backend_wmg is mixed I/O and compute bound

  wait_for_steady_state        = var.wait_for_steady_state
  dd_key_secret_arn            = "arn:aws:secretsmanager:us-west-2:231426846575:secret:dd_api_key-tvi1Ey"
}
