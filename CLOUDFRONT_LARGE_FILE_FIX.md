# CloudFront Large File Caching Fix - Implementation Summary

## Problem Statement

CloudFront has file size caching limits (around 50GB), causing issues when attempting to cache and serve large dataset files. This leads to timeout errors and poor user experience when downloading files larger than 50GB.

## Solution Overview

Configure CloudFront to respect S3 object metadata Cache-Control headers, and automatically set `Cache-Control: no-cache` on files larger than 50GB to bypass CloudFront caching for these large files.

## Changes Implemented

### 1. Infrastructure Changes (Terraform)

**File:** `single-cell-infra/terraform/modules/corpora/backend/s3.tf`

Updated three CloudFront distributions to forward and respect Cache-Control headers:

#### Changes Made:

- Added `headers = ["Cache-Control"]` to `forwarded_values` block
- Added comments explaining the caching behavior
- Applied to all three distributions:
  - `s3_distribution` (datasets.cellxgene.${domain})
  - `s3_wmg_distribution` (ge-data.cellxgene.${domain})
  - `s3_cellguide_distribution` (cellguide.cellxgene.${domain})

#### How It Works:

With `min_ttl = 0` and Cache-Control header forwarding:

- Files WITHOUT Cache-Control metadata: Cached for 1 year (default_ttl)
- Files WITH `Cache-Control: no-cache`: Bypass CloudFront cache, served directly from S3

### 2. Application Changes (Python)

**File:** `backend/layers/thirdparty/s3_provider.py`

#### Updated Methods:

##### `upload_file()` Method

- Automatically checks file size before upload
- For files > 50GB: Sets `CacheControl: "no-cache, no-store, must-revalidate"` in ExtraArgs
- Logs when large files are detected with size information
- Preserves all other upload parameters (ACL, etc.)

##### `upload_directory()` Method

- Scans directory for files > 50GB before upload
- Large files uploaded individually via boto3 with Cache-Control metadata
- Small files uploaded normally via AWS CLI
- Uses exclude patterns to avoid double-uploading large files

### 3. Retroactive Fix Script

**File:** `backend/scripts/update_large_file_cache_control.py`

A comprehensive script to update existing large files in S3:

**Features:**

- Scans specified S3 buckets for objects > 50GB
- Updates metadata without re-uploading files (copy-to-self with new metadata)
- Supports dry-run mode for safe testing
- Per-bucket or stage-wide processing
- Detailed logging and summary statistics
- Preserves all existing object metadata and attributes

**Usage:**

```bash
# Dry run
python update_large_file_cache_control.py --deployment-stage dev --dry-run

# Production update
python update_large_file_cache_control.py --deployment-stage prod --profile single-cell-prod
```

### 4. Documentation

**File:** `backend/scripts/CACHE_CONTROL_LARGE_FILES.md`

Comprehensive documentation covering:

- Architecture overview
- Usage instructions for new and existing files
- Verification procedures
- Troubleshooting guide
- Terraform deployment instructions

## Deployment Steps

### Step 1: Apply Infrastructure Changes

```bash
cd single-cell-infra/terraform/envs/<env>/backend

# Review changes
terraform plan

# Apply (will take 15-30 minutes for CloudFront update)
terraform apply
```

**Environments to update:**

- [ ] dev
- [ ] staging
- [ ] prod

### Step 2: Deploy Application Code

The S3Provider changes will automatically apply to new file uploads once the code is deployed.

```bash
# Standard deployment process for data portal backend
```

### Step 3: Update Existing Large Files

After CloudFront and application updates are deployed:

```bash
cd backend/scripts

# 1. Dry run first (safe to test)
python update_large_file_cache_control.py --deployment-stage dev --dry-run

# 2. Apply to dev
python update_large_file_cache_control.py --deployment-stage dev

# 3. Test in dev, then proceed to staging and prod
python update_large_file_cache_control.py --deployment-stage staging --profile single-cell-dev
python update_large_file_cache_control.py --deployment-stage prod --profile single-cell-prod
```

### Step 4: Verification

```bash
# 1. Check S3 metadata
aws s3api head-object \
  --bucket dataset-assets-public-prod \
  --key <large-file-key> \
  --query 'CacheControl'

# Expected: "no-cache, no-store, must-revalidate"

# 2. Test CloudFront behavior
curl -I https://datasets.cellxgene.cziscience.com/<large-file-key>

# Look for: Cache-Control header and X-Cache: Miss
```

## Benefits

1. **Fixes Caching Issues**: Large files (> 50GB) bypass CloudFront, avoiding size limit errors
2. **Automatic**: New uploads automatically set correct metadata
3. **Transparent**: No changes needed to upload code outside of s3_provider
4. **Configurable**: 50GB threshold can be adjusted if needed
5. **Backward Compatible**: Small files still cached efficiently (1 year TTL)

## Performance Impact

### Large Files (> 50GB)

- **Before**: CloudFront cache timeouts, failed downloads
- **After**: Direct S3 downloads, bypass cache, reliable delivery
- **Trade-off**: Slightly higher latency (no edge caching), but downloads succeed

### Small Files (≤ 50GB)

- **No change**: Still cached at CloudFront edge for fast delivery

## Monitoring

Monitor these metrics after deployment:

1. **CloudFront Cache Hit Ratio**: Should remain high for small files
2. **Download Success Rate**: Should improve for large files
3. **CloudFront 504 Errors**: Should decrease or eliminate
4. **S3 Request Rates**: May increase slightly (direct downloads for large files)

## Rollback Plan

If issues occur:

1. **Terraform**: Revert CloudFront configuration

   ```bash
   git revert <commit>
   terraform apply
   ```

2. **Application**: Revert s3_provider.py changes

   ```bash
   git revert <commit>
   # Redeploy application
   ```

3. **S3 Metadata**: Re-run script to remove Cache-Control
   ```python
   # Modify script to remove Cache-Control instead of adding it
   ```

## Testing Checklist

- [ ] Terraform plan shows only expected CloudFront changes
- [ ] Terraform apply completes successfully
- [ ] Application deployment successful
- [ ] Unit tests pass (if any added)
- [ ] Manual test: Upload file < 50GB → no Cache-Control set
- [ ] Manual test: Upload file > 50GB → Cache-Control set
- [ ] Script dry-run identifies correct files
- [ ] Script successfully updates metadata in dev
- [ ] CloudFront serves small files with cache hits
- [ ] CloudFront serves large files without caching
- [ ] No 504 errors for large file downloads

## Files Changed

### Infrastructure

- `single-cell-infra/terraform/modules/corpora/backend/s3.tf`

### Application

- `backend/layers/thirdparty/s3_provider.py`

### Scripts & Documentation

- `backend/scripts/update_large_file_cache_control.py` (new)
- `backend/scripts/CACHE_CONTROL_LARGE_FILES.md` (new)
- `CLOUDFRONT_LARGE_FILE_FIX.md` (this file, new)

## References

- AWS CloudFront: [Caching Based on Headers](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/header-caching.html)
- AWS S3: [Object Metadata](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingMetadata.html)
- HTTP Cache-Control: [MDN Web Docs](https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Cache-Control)

## Support

For questions or issues:

1. Check `backend/scripts/CACHE_CONTROL_LARGE_FILES.md` for detailed troubleshooting
2. Review CloudFront logs in `<bucket>-logs/cloudfront/`
3. Check S3 access logs for direct downloads
4. Contact infrastructure team for CloudFront distribution issues
