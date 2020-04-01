locals {
  tags = {
    project   = var.project
    env       = var.env
    service   = var.service
    owner     = var.owner
    managedBy = "terraform"
  }

  domain       = replace(data.aws_route53_zone.zone.name, "/\\.$/", "")
  website_fqdn = "${var.subdomain}.${local.domain}"
  bucket_name  = var.bucket_name != "" ? var.bucket_name : local.website_fqdn

  aliases = [
    "${local.website_fqdn}",
    "www.${local.website_fqdn}",
  ]
}

data "aws_route53_zone" "zone" {
  zone_id = var.aws_route53_zone_id
}

resource "aws_cloudfront_origin_access_identity" "origin_access_identity" {}

data "aws_iam_policy_document" "bucket_policy" {
  statement {
    sid       = "getObjects"
    actions   = ["s3:GetObject"]
    resources = ["arn:aws:s3:::${local.bucket_name}/*"]
    principals {
      type        = "AWS"
      identifiers = ["${aws_cloudfront_origin_access_identity.origin_access_identity.iam_arn}"]
    }
//    principals {
//      type = "*"
//      identifiers = [
//      "*"]
//    }
  }

  statement {
    sid       = "listBucket"
    actions   = ["s3:ListBucket"]
    resources = ["arn:aws:s3:::${local.bucket_name}"]
    principals {
      type        = "AWS"
      identifiers = ["${aws_cloudfront_origin_access_identity.origin_access_identity.iam_arn}"]
    }
//    principals {
//      type = "*"
//      identifiers = [
//      "*"]
//    }
  }
}

resource "aws_s3_bucket" "bucket" {
  bucket = local.bucket_name
  acl    = "private"
  policy = data.aws_iam_policy_document.bucket_policy.json

  // Cloudfront needs this to compress assets
  // https://stackoverflow.com/questions/35590622/cloudfront-with-s3-website-as-origin-is-not-serving-gzipped-files
  // Content-Type is also needed to allow CORS json requests
  cors_rule {
    allowed_headers = ["Authorization", "Content-Length", "Content-Type"]
    allowed_methods = ["GET"]
    allowed_origins = ["*"]
    max_age_seconds = 3000
  }

  server_side_encryption_configuration {
    rule {
      apply_server_side_encryption_by_default {
        sse_algorithm = "AES256"
      }
    }
  }

  website {
    index_document = "index.html"
    error_document = "index.html"
  }
  tags = local.tags
}

resource "aws_s3_bucket_public_access_block" "bucket" {
  bucket = aws_s3_bucket.bucket.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_cloudfront_distribution" "s3_distribution" {
  origin {
//    domain_name = aws_s3_bucket.bucket.website_domain
    domain_name = aws_s3_bucket.bucket.bucket_domain_name
    origin_id   = local.website_fqdn

    s3_origin_config {
      origin_access_identity = aws_cloudfront_origin_access_identity.origin_access_identity.cloudfront_access_identity_path
    }
  }

  enabled             = true
  is_ipv6_enabled     = true
  default_root_object = var.index_document_path

  aliases = concat(var.aliases, local.aliases)

  default_cache_behavior {
    allowed_methods = ["GET", "HEAD", "OPTIONS"]
    cached_methods  = ["GET", "HEAD"]

    target_origin_id = local.website_fqdn

    forwarded_values {
      query_string = false

      cookies {
        forward = "none"
      }
    }

    viewer_protocol_policy = "redirect-to-https"
    min_ttl                = 0
    default_ttl            = 3600
    max_ttl                = 86400
    compress               = true

    lambda_function_association {
        event_type= "viewer-request"
        lambda_arn="arn:aws:lambda:us-east-1:699936264352:function:single-cell-site-index-forwarding-test:1"
        include_body = false
    }
  }

  ordered_cache_behavior {
    path_pattern    = var.path_pattern
    allowed_methods = ["GET", "HEAD", "OPTIONS"]
    cached_methods  = ["GET", "HEAD"]

    target_origin_id = local.website_fqdn

    forwarded_values {
      query_string = true
      headers      = ["Origin"]

      cookies {
        forward = "none"
      }
    }

    viewer_protocol_policy = "redirect-to-https"
    min_ttl                = 0
    default_ttl            = 3600
    max_ttl                = 86400
    compress               = true
  }

  restrictions {
    geo_restriction {
      restriction_type = "none"
    }
  }

  price_class = var.cloudfront_price_class

  viewer_certificate {
    acm_certificate_arn      = var.aws_acm_cert_arn
    ssl_support_method       = "sni-only"
    minimum_protocol_version = var.minimum_tls_version
  }

  # This is error handling logic for single page applications
  custom_error_response {
    error_code         = 403
    response_code      = 200
    response_page_path = "/${var.index_document_path}"
  }

  custom_error_response {
    error_code         = 404
    response_code      = 200
    response_page_path = "/${var.index_document_path}"
  }

  custom_error_response {
    error_code         = 503
    response_code      = 200
    response_page_path = "/${var.index_document_path}"
  }

  tags = local.tags
}

# for ipv4
resource "aws_route53_record" "ipv4-record" {
  zone_id = var.aws_route53_zone_id
  name    = local.website_fqdn
  type    = "A"

  alias {
    name                   = aws_cloudfront_distribution.s3_distribution.domain_name
    zone_id                = aws_cloudfront_distribution.s3_distribution.hosted_zone_id
    evaluate_target_health = true
  }
}

# for ipv6
resource "aws_route53_record" "ipv6-record" {
  zone_id = var.aws_route53_zone_id
  name    = local.website_fqdn
  type    = "AAAA"

  alias {
    name                   = aws_cloudfront_distribution.s3_distribution.domain_name
    zone_id                = aws_cloudfront_distribution.s3_distribution.hosted_zone_id
    evaluate_target_health = true
  }
}

# for ipv4
resource "aws_route53_record" "www-ipv4-record" {
  zone_id = var.aws_route53_zone_id
  name    = "www.${local.website_fqdn}"
  type    = "A"

  alias {
    name                   = aws_cloudfront_distribution.s3_distribution.domain_name
    zone_id                = aws_cloudfront_distribution.s3_distribution.hosted_zone_id
    evaluate_target_health = true
  }
}

# for ipv6
resource "aws_route53_record" "www-ipv6-record" {
  zone_id = var.aws_route53_zone_id
  name    = "www.${local.website_fqdn}"
  type    = "AAAA"

  alias {
    name                   = aws_cloudfront_distribution.s3_distribution.domain_name
    zone_id                = aws_cloudfront_distribution.s3_distribution.hosted_zone_id
    evaluate_target_health = true
  }
}
