resource "aws_secretsmanager_secret" "browser_database_secrets" {
  name = "dcp/backend/browser/${var.deployment_stage}/database"
}

resource "aws_secretsmanager_secret_version" "browser_database_secrets" {
  secret_id = aws_secretsmanager_secret.browser_database_secrets.id
  secret_string = <<SECRETS_JSON
{
  "database_uri": "mysql+pymysql://${aws_rds_cluster.browser.master_username}:${aws_rds_cluster.browser.master_password}@${aws_rds_cluster.browser.endpoint}/${aws_rds_cluster.browser.database_name}"
}
SECRETS_JSON
}
