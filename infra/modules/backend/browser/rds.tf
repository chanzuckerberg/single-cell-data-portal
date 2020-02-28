resource "aws_rds_cluster_instance" "cluster_instances" {
  count                        = var.db_instance_count
  identifier                   = "browser-cluster-${var.deployment_stage}-${count.index}"
  cluster_identifier           = aws_rds_cluster.browser.id
  instance_class               = "db.r4.large"
  publicly_accessible          = "true"
  engine                       = "aurora-mysql"
  engine_version               = "5.7.mysql_aurora.2.07.1"
  auto_minor_version_upgrade   = "true"
  performance_insights_enabled = "true"
  preferred_maintenance_window = var.preferred_maintenance_window
}

resource "aws_rds_cluster" "browser" {
  cluster_identifier              = "dcp-browser-${var.deployment_stage}"
  engine                          = "aurora-mysql"
  engine_version                  = "5.7.mysql_aurora.2.07.1"
  availability_zones              = ["us-east-1a", "us-east-1c", "us-east-1d"]
  database_name                   = "browser_${var.deployment_stage}"
  master_username                 = var.db_username
  master_password                 = var.db_password
  backup_retention_period         = 7
  port                            = 3306
  preferred_backup_window         = "07:27-07:57"
  preferred_maintenance_window    = var.preferred_maintenance_window
  storage_encrypted               = "true"
  vpc_security_group_ids          = [aws_security_group.rds_mysql.id]
  db_cluster_parameter_group_name = "default.aurora-mysql5.7"
  final_snapshot_identifier       = "dcp-browser-${var.deployment_stage}"
}

resource "aws_default_vpc" "default" {
  tags = {
    Name = "Default VPC"
  }
}

resource "aws_security_group" "rds_mysql" {
  name        = "dcp-browser-${var.deployment_stage}-rds-mysql-sg"
  description = "DCP browser rds security group"
  vpc_id      = aws_default_vpc.default.id

  ingress {
    from_port   = 3306
    to_port     = 3306
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
    description = "all"
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}
