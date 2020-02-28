resource "aws_secretsmanager_secret" "ledger_database_secrets" {
  name = "dcp/backend/ledger/${var.deployment_stage}/database"
}

resource "aws_secretsmanager_secret_version" "ledger_database_secrets" {
  secret_id     = aws_secretsmanager_secret.ledger_database_secrets.id
  secret_string = <<SECRETS_JSON
{
  "database_uri": "postgresql://${aws_rds_cluster.ledger.master_username}:${aws_rds_cluster.ledger.master_password}@${aws_rds_cluster.ledger.endpoint}/${aws_rds_cluster.ledger.database_name}"
}
SECRETS_JSON
}