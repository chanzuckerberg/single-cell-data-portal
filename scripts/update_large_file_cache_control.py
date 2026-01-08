#!/usr/bin/env python3
"""
Script to retroactively update Cache-Control metadata on large S3 objects (exceeding CloudFront max cache size)
to prevent CloudFront caching issues.

This script:
1. Scans specified S3 buckets for objects larger than CloudFront max cache size
2. Updates their metadata to set Cache-Control: no-cache, no-store, must-revalidate
3. Logs all updates for audit purposes

The CloudFront max cache size can be configured via the CLOUDFRONT_MAX_CACHE_SIZE environment variable
(e.g., "50GB", "100MB", "1TB"). Defaults to 50GB if not set.

Usage:
    python update_large_file_cache_control.py --deployment-stage <dev|staging|prod> [--dry-run] [--bucket <bucket-name>]

Examples:
    # Dry run to see what would be updated in dev
    python update_large_file_cache_control.py --deployment-stage dev --dry-run
    
    # Update all buckets in staging
    python update_large_file_cache_control.py --deployment-stage staging
    
    # Update only dataset-assets-public bucket in prod
    python update_large_file_cache_control.py --deployment-stage prod --bucket dataset-assets-public-prod
"""

import argparse
import logging
import sys
from typing import List, Tuple

import boto3
from botocore.exceptions import ClientError

from backend.layers.thirdparty.s3_provider import CLOUDFRONT_MAX_CACHE_SIZE

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# Constants
CACHE_CONTROL_VALUE = "no-cache, no-store, must-revalidate"


def get_buckets_for_stage(deployment_stage: str) -> List[str]:
    """
    Returns the list of buckets to scan for a given deployment stage.

    Args:
        deployment_stage: The deployment stage (dev, staging, prod)

    Returns:
        List of bucket names
    """
    return [
        f"dataset-assets-public-{deployment_stage}",
        f"gene-expression-assets-public-{deployment_stage}",
        f"cellguide-data-public-{deployment_stage}",
    ]


def find_large_objects(bucket_name: str, s3_client) -> List[Tuple[str, int]]:
    """
    Finds all objects in a bucket that are larger than CloudFront max cache size.

    Args:
        bucket_name: The name of the S3 bucket
        s3_client: Boto3 S3 client

    Returns:
        List of tuples (object_key, size_in_bytes)
    """
    large_objects = []

    logger.info(f"Scanning bucket: {bucket_name}")

    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        page_iterator = paginator.paginate(Bucket=bucket_name)

        total_objects = 0
        for page in page_iterator:
            if "Contents" not in page:
                continue

            for obj in page["Contents"]:
                total_objects += 1
                if obj["Size"] > CLOUDFRONT_MAX_CACHE_SIZE:
                    large_objects.append((obj["Key"], obj["Size"]))
                    size_gb = obj["Size"] / (1024 * 1024 * 1024)
                    logger.info(f"Found large object: {obj['Key']} ({size_gb:.2f} GB)")

        logger.info(f"Scanned {total_objects} objects in {bucket_name}, found {len(large_objects)} large objects")

    except ClientError as e:
        logger.error(f"Error scanning bucket {bucket_name}: {e}")
        raise

    return large_objects


def update_object_metadata(bucket_name: str, object_key: str, s3_client, dry_run: bool = False) -> bool:
    """
    Updates the Cache-Control metadata for a single S3 object.

    Args:
        bucket_name: The name of the S3 bucket
        object_key: The key of the object to update
        s3_client: Boto3 S3 client
        dry_run: If True, only log what would be done without making changes

    Returns:
        True if successful, False otherwise
    """
    try:
        # First, get the current object metadata
        head_response = s3_client.head_object(Bucket=bucket_name, Key=object_key)

        # Check if Cache-Control is already set correctly
        current_cache_control = head_response.get("CacheControl", "")
        if current_cache_control == CACHE_CONTROL_VALUE:
            logger.info(f"Object {object_key} already has correct Cache-Control metadata, skipping")
            return True

        if dry_run:
            logger.info(f"[DRY RUN] Would update Cache-Control for {object_key}")
            logger.info(f"[DRY RUN] Current: '{current_cache_control}' -> New: '{CACHE_CONTROL_VALUE}'")
            return True

        # Copy the object to itself with new metadata
        # This is the standard way to update metadata in S3
        copy_source = {"Bucket": bucket_name, "Key": object_key}

        # Preserve existing metadata while adding/updating Cache-Control
        metadata_directive = "REPLACE"
        metadata = head_response.get("Metadata", {})

        # Build ExtraArgs for copy operation
        extra_args = {
            "CacheControl": CACHE_CONTROL_VALUE,
            "Metadata": metadata,
            "MetadataDirective": metadata_directive,
        }

        # Preserve ACL if present
        if "ACL" in head_response:
            extra_args["ACL"] = "bucket-owner-full-control"

        # Preserve content type
        if "ContentType" in head_response:
            extra_args["ContentType"] = head_response["ContentType"]

        # Preserve content encoding
        if "ContentEncoding" in head_response:
            extra_args["ContentEncoding"] = head_response["ContentEncoding"]

        # Perform the copy operation
        s3_client.copy_object(
            Bucket=bucket_name,
            Key=object_key,
            CopySource=copy_source,
            **extra_args,
        )

        logger.info(f"Successfully updated Cache-Control for {object_key}")
        return True

    except ClientError as e:
        logger.error(f"Error updating object {object_key}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Update Cache-Control metadata for large S3 objects (exceeding CloudFront max cache size)"
    )
    parser.add_argument(
        "--deployment-stage",
        required=True,
        choices=["dev", "staging", "prod"],
        help="Deployment stage (dev, staging, or prod)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run without making actual changes",
    )
    parser.add_argument(
        "--bucket",
        help="Specific bucket to process (optional, default: all buckets for stage)",
    )
    parser.add_argument(
        "--profile",
        help="AWS profile to use (optional)",
    )

    args = parser.parse_args()

    # Initialize S3 client
    session_args = {}
    if args.profile:
        session_args["profile_name"] = args.profile

    session = boto3.Session(**session_args)
    s3_client = session.client("s3")

    # Determine which buckets to process
    buckets = [args.bucket] if args.bucket else get_buckets_for_stage(args.deployment_stage)

    logger.info(f"{'[DRY RUN] ' if args.dry_run else ''}Starting Cache-Control metadata update")
    logger.info(f"Deployment stage: {args.deployment_stage}")
    logger.info(f"Buckets to process: {buckets}")

    # Process each bucket
    total_large_objects = 0
    total_updated = 0
    total_failed = 0

    for bucket_name in buckets:
        logger.info(f"\n{'='*80}")
        logger.info(f"Processing bucket: {bucket_name}")
        logger.info(f"{'='*80}\n")

        try:
            # Find large objects
            large_objects = find_large_objects(bucket_name, s3_client)
            total_large_objects += len(large_objects)

            # Update metadata for each large object
            for object_key, size in large_objects:
                size_gb = size / (1024 * 1024 * 1024)
                logger.info(f"\nProcessing: {object_key} ({size_gb:.2f} GB)")

                success = update_object_metadata(bucket_name, object_key, s3_client, args.dry_run)
                if success:
                    total_updated += 1
                else:
                    total_failed += 1

        except Exception as e:
            logger.error(f"Error processing bucket {bucket_name}: {e}")
            continue

    # Print summary
    logger.info(f"\n{'='*80}")
    logger.info("Summary")
    logger.info(f"{'='*80}")
    logger.info(f"Total large objects found: {total_large_objects}")
    logger.info(f"Successfully updated: {total_updated}")
    logger.info(f"Failed: {total_failed}")

    if args.dry_run:
        logger.info("\nThis was a DRY RUN. No changes were made to S3.")
        logger.info("Run without --dry-run to apply the changes.")

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
