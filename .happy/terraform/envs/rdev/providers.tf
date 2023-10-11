provider aws {
  version = "~> 4.23.0"
  region  = "us-west-2"
  assume_role {
    role_arn = "arn:aws:iam::${var.aws_account_id}:role/${var.aws_role}"
  }
  allowed_account_ids = [var.aws_account_id]
}
