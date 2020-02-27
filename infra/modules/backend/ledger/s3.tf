variable "bucket_prefix" {
  type = string
  default = "dcp-ledger-bucket"
}

resource "aws_s3_bucket" "ledger_s3_bucket" {
  bucket = "${var.bucket_prefix}-${var.deployment_stage}"
  acl    = "private"

  logging {
    target_bucket = var.log_bucket
    target_prefix = "log/${var.bucket_prefix}-${var.deployment_stage}"
  }

  // All s3 buckets should have this enabled
  server_side_encryption_configuration {
    rule {
      apply_server_side_encryption_by_default {
        sse_algorithm     = "AES256"
      }
    }
  }
}
