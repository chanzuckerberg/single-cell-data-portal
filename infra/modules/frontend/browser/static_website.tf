resource "aws_s3_bucket" "gatsby_static_bucket" {

  bucket = "dcp-site-deployment-${var.deployment_stage}"

  website {
    index_document = "index.html"
    error_document = "404.html"
  }
}

resource "aws_s3_bucket_public_access_block" "gatsby_static_bucket_publicaccess" {
  bucket = "${aws_s3_bucket.gatsby_static_bucket.id}"

  block_public_acls       = true
  block_public_policy     = false
  ignore_public_acls      = true
  restrict_public_buckets = false
}

data "aws_iam_policy_document" "gatsby_static_bucket_policy_document" {
  statement {
    actions = [
      "s3:GetObject"
    ]

    resources = [
      "${aws_s3_bucket.gatsby_static_bucket.arn}/*"
    ]

    principals {
      type = "*"
      identifiers = [
      "*"]
    }

    effect = "Allow"
  }
}

resource "aws_s3_bucket_policy" "gatsby_static_bucket_policy" {
  bucket = "${aws_s3_bucket.gatsby_static_bucket.id}"
  policy = "${data.aws_iam_policy_document.gatsby_static_bucket_policy_document.json}"
}
