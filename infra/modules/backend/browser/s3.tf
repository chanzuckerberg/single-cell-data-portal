resource "aws_s3_bucket" "browser_s3_bucket" {
  bucket = "dcp-browser-bucket-${var.deployment_stage}"
  acl    = "private"
}
