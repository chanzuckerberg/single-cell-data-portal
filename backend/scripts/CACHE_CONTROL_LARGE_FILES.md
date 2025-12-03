# Cache-Control Configuration for Large S3 Files

## Overview

CloudFront has file size caching limits (around 50GB for standard configurations). To work around this limitation, we configure CloudFront to respect S3 object metadata Cache-Control headers, and automatically set `Cache-Control: no-cache` on files larger than 50GB.

## Architecture

### 1. CloudFront Configuration

The CloudFront distributions are configured to:

- Forward the `Cache-Control` header from S3 origins
- Use `min_ttl=0` which allows CloudFront to respect origin Cache-Control headers
- Apply a default TTL of 1 year for files without Cache-Control metadata

**Affected CloudFront distributions:**

- `s3_distribution` (datasets.cellxgene.${domain})
- `s3_wmg_distribution` (ge-data.cellxgene.${domain})
- `s3_cellguide_distribution` (cellguide.cellxgene.${domain})

**Terraform location:** `single-cell-infra/terraform/modules/corpora/backend/s3.tf`

### 2. Automatic Metadata Setting on Upload

The S3Provider class automatically sets Cache-Control metadata during file uploads:

- **Location:** `backend/layers/thirdparty/s3_provider.py`
- **Methods updated:**
  - `upload_file()` - Checks file size and sets Cache-Control for files > 50GB
  - `upload_directory()` - Identifies large files and uploads them with proper metadata

**Behavior:**

- Files ≤ 50GB: Uploaded normally, cached by CloudFront (default 1 year TTL)
- Files > 50GB: Uploaded with `Cache-Control: no-cache, no-store, must-revalidate`

### 3. Cache-Control Value

For large files, we use:

```
Cache-Control: no-cache, no-store, must-revalidate
```

This ensures:

- `no-cache`: CloudFront must revalidate with origin before serving cached copy
- `no-store`: Content should not be stored in any cache
- `must-revalidate`: Forces caches to check with origin before serving stale content

## Usage

### For New Uploads

No action required! The S3Provider automatically handles Cache-Control metadata:

```python
from backend.layers.thirdparty.s3_provider import S3Provider

s3_provider = S3Provider()

# This will automatically set Cache-Control if file is > 50GB
s3_provider.upload_file(
    src_file="/path/to/large/file.h5ad",
    bucket_name="dataset-assets-public-dev",
    dst_file="dataset-id/file.h5ad",
    extra_args={"ACL": "bucket-owner-full-control"}
)
```

### For Existing Large Files

Use the provided script to retroactively update Cache-Control metadata on existing large files:

#### Dry Run (Recommended First)

```bash
cd backend/scripts

# Check what would be updated in dev environment
python update_large_file_cache_control.py --deployment-stage dev --dry-run

# Check staging
python update_large_file_cache_control.py --deployment-stage staging --dry-run --profile single-cell-dev
```

#### Update Production

```bash
# Update all buckets in production
python update_large_file_cache_control.py --deployment-stage prod --profile single-cell-prod

# Update only specific bucket
python update_large_file_cache_control.py \
  --deployment-stage prod \
  --bucket dataset-assets-public-prod \
  --profile single-cell-prod
```

#### Script Options

- `--deployment-stage`: Required. Environment to process (dev, staging, prod)
- `--dry-run`: Optional. Preview changes without modifying anything
- `--bucket`: Optional. Process only a specific bucket (default: all stage buckets)
- `--profile`: Optional. AWS profile to use for credentials

## Verification

### 1. Check S3 Object Metadata

```bash
# Using AWS CLI
aws s3api head-object \
  --bucket dataset-assets-public-dev \
  --key <object-key> \
  --query 'CacheControl'

# Expected output for large files:
# "no-cache, no-store, must-revalidate"
```

### 2. Check CloudFront Behavior

```bash
# Check headers returned by CloudFront
curl -I https://datasets.cellxgene.cziscience.com/<object-key>

# Look for:
# Cache-Control: no-cache, no-store, must-revalidate
# X-Cache: Miss from cloudfront (for uncached files)
```

### 3. Monitor CloudFront Logs

CloudFront logs are stored in:

- `dataset-assets-public-${env}-logs/cloudfront/`

Look for:

- Hit/Miss ratios for large files
- Any 504 Gateway Timeout errors (indicating caching issues)

## Terraform Changes

To apply the CloudFront configuration changes:

```bash
cd terraform/envs/<env>/backend

# Plan the changes
terraform plan

# Apply the changes
terraform apply

# Note: CloudFront distribution updates can take 15-30 minutes to deploy
```

## Troubleshooting

### Issue: Large files still timing out

**Possible causes:**

1. CloudFront distribution hasn't been updated yet
2. Object metadata not set correctly
3. CloudFront cache contains old cached versions

**Solutions:**

```bash
# 1. Verify CloudFront configuration
aws cloudfront get-distribution-config --id <distribution-id>

# 2. Check object metadata
aws s3api head-object --bucket <bucket> --key <key>

# 3. Invalidate CloudFront cache for specific files
aws cloudfront create-invalidation \
  --distribution-id <distribution-id> \
  --paths "/<object-key>"
```

### Issue: Script fails with access denied

**Solution:**
Ensure your AWS credentials have the following permissions:

- `s3:ListBucket`
- `s3:GetObject`
- `s3:GetObjectMetadata`
- `s3:PutObject` (for metadata updates)

### Issue: Metadata update is slow for many files

The script processes files sequentially to avoid overwhelming S3 API limits. For large-scale updates:

1. Run the script on a powerful EC2 instance in the same region
2. Consider processing one bucket at a time
3. Use `--bucket` option to process specific buckets

## Implementation Timeline

1. ✅ Update CloudFront distributions (Terraform)
2. ✅ Update S3Provider upload methods
3. ⏳ Apply Terraform changes to environments
4. ⏳ Run retroactive script on existing files
5. ⏳ Verify no caching issues with large files

## References

- [CloudFront Caching Based on Request Headers](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/header-caching.html)
- [S3 Object Metadata](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingMetadata.html)
- [CloudFront Cache Behavior Settings](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/distribution-web-values-specify.html)
