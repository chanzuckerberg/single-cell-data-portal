resource "aws_s3_bucket" "s3_log_bucket" {
  bucket = "dcp-log-bucket-${var.deployment_stage}"
  acl    = "log-delivery-write"

  lifecycle_rule {
    id      = "log"
    enabled = true

    prefix = "log/"

    tags = {
      "rule"      = "log"
      "autoclean" = "true"
    }

    expiration {
      days = var.log_retention_time
    }
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

output "log_bucket" {
  value = aws_s3_bucket.s3_log_bucket.id
}
