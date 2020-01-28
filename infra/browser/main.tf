#make bucket
resource "aws_s3_bucket" "browser" {
  bucket_prefix = "dcp-browser-${var.DCP_DEPLOYMENT_STAGE}"
  acl    = "public-read"
  force_destroy = true
  server_side_encryption_configuration {
    rule {
      apply_server_side_encryption_by_default {
        sse_algorithm     = "AES256"
      }
    }
  }
  tags = {
    Environment = "${var.DCP_DEPLOYMENT_STAGE}"
  }
}

#TODO make policy

#TODO make api gateway


