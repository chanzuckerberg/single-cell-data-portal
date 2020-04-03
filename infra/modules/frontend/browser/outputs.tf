output "bucket_name" {
  value = local.bucket_name
}

output "bucket_arn" {
  value = aws_s3_bucket.bucket.arn
}

output "cloudfront_arn" {
  value = aws_cloudfront_distribution.s3_distribution.arn
}

output "cloudfront_domain_name" {
  value = aws_cloudfront_distribution.s3_distribution.domain_name
}

output "cloudfront_hosted_zone_id" {
  value = aws_cloudfront_distribution.s3_distribution.hosted_zone_id
}
