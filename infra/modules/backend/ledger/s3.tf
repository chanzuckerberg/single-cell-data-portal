resource "aws_s3_bucket" "ledger_s3_bucket" {
  bucket = "dcp-ledger-bucket-${var.deployment_stage}"
  acl    = "private"
}