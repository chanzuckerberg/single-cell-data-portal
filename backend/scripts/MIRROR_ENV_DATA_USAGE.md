# Mirror Environment Data - Usage Guide

This guide explains how to mirror RDS database and S3 data between environments.

## Usage

The `mirror_env_data` recipe ([Makefile#L137](../Makefile#L137)) automatically manages SSH tunnels for you!

### Basic Command

```bash
cd backend
make mirror_env_data SRC_ENV=prod DEST_ENV=staging
```

### Common Examples

> **⚠️ Always test with `DRY_RUN=1` first!** See [Dry-Run Mode](#dry-run-mode) for details.

```bash
# Test first (recommended) - shows what would be changed without modifying data
make mirror_env_data SRC_ENV=prod DEST_ENV=staging DRY_RUN=1

# Mirror prod to staging (most common use case)
make mirror_env_data SRC_ENV=prod DEST_ENV=staging

# Mirror prod to dev
make mirror_env_data SRC_ENV=prod DEST_ENV=dev

# Mirror with WMG cube
make mirror_env_data SRC_ENV=prod DEST_ENV=staging WMG_CUBE=1

# Mirror with CellGuide data to rdev
make mirror_env_data SRC_ENV=dev DEST_ENV=rdev STACK=my-stack CELLGUIDE=1

# Mirror specific collections to rdev (with sample data)
make mirror_env_data SRC_ENV=dev DEST_ENV=rdev STACK=my-stack COLLECTIONS=uuid1,uuid2

# Mirror specific collections to rdev (with real data)
make mirror_env_data SRC_ENV=dev DEST_ENV=rdev STACK=my-stack COLLECTIONS=uuid1,uuid2 DATA=1
```

### Parameters

| Parameter     | Description                                                            | Default             | Script Reference                                     |
| ------------- | ---------------------------------------------------------------------- | ------------------- | ---------------------------------------------------- |
| `SRC_ENV`     | Source environment (`prod`, `staging`, `dev`)                          | `prod`              | [Makefile:133](../Makefile#L133)                     |
| `DEST_ENV`    | Destination environment (`dev`, `staging`, `rdev`)                     | Required            | [set_src_dest_envs.sh#L5](set_src_dest_envs.sh#L5)   |
| `DRY_RUN`     | Test mode - show what would change without modifying data (set to `1`) | -                   | [mirror_env_data.sh#L19](mirror_env_data.sh#L19)     |
| `WMG_CUBE`    | Copy WMG cube (set to `1` to enable)                                   | -                   | [mirror_s3_data.sh#L36](mirror_s3_data.sh#L36)       |
| `CELLGUIDE`   | Copy CellGuide data (set to `1` to enable, `rdev` only)                | -                   | [mirror_s3_data.sh#L51](mirror_s3_data.sh#L51)       |
| `STACK`       | Stack name for rdev environment                                        | Required for `rdev` | [set_src_dest_envs.sh#L22](set_src_dest_envs.sh#L22) |
| `COLLECTIONS` | Comma-separated Collection UUIDs (`rdev` only)                         | -                   | [mirror_s3_data.sh#L62](mirror_s3_data.sh#L62)       |
| `DATA`        | Use real datafiles instead of sample data (`rdev` only, set to `1`)    | -                   | [mirror_s3_data.sh#L94](mirror_s3_data.sh#L94)       |
| `SRC_PORT`    | Override source tunnel port                                            | `5432`              | [Makefile:134](../Makefile#L134)                     |
| `DEST_PORT`   | Override destination tunnel port                                       | `5433`              | [Makefile:135](../Makefile#L135)                     |

## How It Works

The mirroring process runs in two phases ([mirror_env_data.sh](mirror_env_data.sh)):

### Phase 1: Database Mirroring ([mirror_rds_data.sh](mirror_rds_data.sh))

1. **Validates environments** ([set_src_dest_envs.sh#L16](set_src_dest_envs.sh#L16)) and prompts for confirmation
2. **Sets up SSH tunnels** using AWS SSM ([manage_db_tunnel.sh#L17](manage_db_tunnel.sh#L17))
3. **Dumps source database** using `pg_dump` via Makefile ([mirror_rds_data.sh#L26](mirror_rds_data.sh#L26))
4. **Backs up destination** (unless `DEST_ENV=rdev`) ([mirror_rds_data.sh#L32](mirror_rds_data.sh#L32))
5. **Restores to destination** using `pg_restore` ([mirror_rds_data.sh#L52](mirror_rds_data.sh#L52))
6. **Updates environment-specific URLs** in database ([mirror_rds_data.sh#L59](mirror_rds_data.sh#L59))
7. **Closes tunnels** ([mirror_env_data.sh#L95](mirror_env_data.sh#L95))

### Phase 2: S3 Mirroring ([mirror_s3_data.sh](mirror_s3_data.sh))

- Syncs S3 buckets using AWS CLI (`aws s3 sync` and `aws s3 cp`)
- Handles WMG cube, CellGuide data, and collection-specific data based on parameters
- **For standard environments** (dev, staging):
  - Uses 16 concurrent parallel jobs for faster transfer ([mirror_s3_data.sh#L131](mirror_s3_data.sh#L131))
  - Syncs three main buckets: `corpora-data-*`, `hosted-cellxgene-*`, and `dataset-assets-public-*`
  - Uses hex-based parallelization (0-f) to distribute workload ([mirror_s3_data.sh#L133](mirror_s3_data.sh#L133))
- **For rdev with COLLECTIONS**:
  - Queries database to get Dataset URIs for specified Collections ([mirror_s3_data.sh#L81](mirror_s3_data.sh#L81))
  - Uses GNU Parallel's `sem` with unlimited jobs based on CPU cores ([mirror_s3_data.sh#L95](mirror_s3_data.sh#L95))
  - Copies either sample data or real data based on `DATA` parameter
- Copies metadata without tags using `--copy-props metadata-directive` ([mirror_s3_data.sh#L25](mirror_s3_data.sh#L25))
- Supports dry-run mode via `--dryrun` flag when `DRY_RUN=1` ([mirror_s3_data.sh#L28](mirror_s3_data.sh#L28))

### Dry-Run Mode

When `DRY_RUN=1` is set ([mirror_env_data.sh#L19](mirror_env_data.sh#L19)):

- Sets up tunnels to validate connectivity
- Skips all database dump/restore operations
- Runs S3 sync/copy operations with `--dryrun` flag to show what would be changed without copying files ([mirror_s3_data.sh#L28](mirror_s3_data.sh#L28))

## Supported Mirroring Paths

Environment validation is handled in [set_src_dest_envs.sh](set_src_dest_envs.sh):

| Source    | Destination | Purpose                                  | Status       | Validation                      |
| --------- | ----------- | ---------------------------------------- | ------------ | ------------------------------- |
| `prod`    | `staging`   | Sync staging with production             | ✅ Supported | [L54](set_src_dest_envs.sh#L54) |
| `prod`    | `dev`       | Sync dev with production                 | ✅ Supported | [L54](set_src_dest_envs.sh#L54) |
| `prod`    | `rdev`      | Sync rdev with production                | ✅ Supported | [L54](set_src_dest_envs.sh#L54) |
| `staging` | `dev`       | Testing/debugging staging issues in dev  | ✅ Supported | [L46](set_src_dest_envs.sh#L46) |
| `staging` | `rdev`      | Testing/debugging staging issues in rdev | ✅ Supported | [L46](set_src_dest_envs.sh#L46) |
| `dev`     | `dev`       | Testing the mirror operation             | ✅ Supported | [L30](set_src_dest_envs.sh#L30) |
| `dev`     | `rdev`      | Testing/debugging dev issues in rdev     | ✅ Supported | [L38](set_src_dest_envs.sh#L38) |
| `rdev`    | anywhere    | Not supported                            | ❌ Blocked   | [L16](set_src_dest_envs.sh#L16) |

## Important Warnings

### ⚠️ This Operation is DESTRUCTIVE

- **Source environment**: NEVER modified (read-only dump operations only - see [mirror_rds_data.sh#L26](mirror_rds_data.sh#L26))
- **Destination environment**: Data will be COMPLETELY REPLACED
- **Backup**: Destination database is backed up before replacement (except for `rdev` - see [mirror_rds_data.sh#L32](mirror_rds_data.sh#L32))
- **Confirmation**: You will be prompted before destructive operations begin ([set_src_dest_envs.sh](set_src_dest_envs.sh))
- **S3 Sync**: Adds objects to destination but doesn't delete existing objects (see [mirror_s3_data.sh#L6](mirror_s3_data.sh#L6))

**Always test with `DRY_RUN=1` first!** See [Dry-Run Mode](#dry-run-mode) above.

### Required Permissions

- **AWS profiles**: `single-cell-prod`, `single-cell-dev` ([\_mirror_utils.sh#L15](_mirror_utils.sh#L15))
  - Note: `staging`, `dev`, and `rdev` all use the `single-cell-dev` profile
- **AWS region**: `us-west-2` ([\_mirror_utils.sh#L7](_mirror_utils.sh#L7))
- **Access to**: RDS databases, S3 buckets, SSM (Systems Manager), EC2 instances, and Secrets Manager
- **Required AWS credentials** for both source and destination accounts
- **For dataset syncs**: Script assumes `sync-datasets-dev` IAM role ([mirror_s3_data.sh#L141](mirror_s3_data.sh#L141))

## Troubleshooting

### Port Already in Use

If you get an error about ports being in use, specify different ports:

```bash
make mirror_env_data SRC_ENV=prod DEST_ENV=staging SRC_PORT=5434 DEST_PORT=5435
```

### Tunnel Connection Issues

If you see `TargetNotConnected` errors, check:

```bash
# Verify EC2 instance is running
aws ec2 describe-instances --region us-west-2 --instance-ids <instance-id>

# Verify SSM connection
aws ssm describe-instance-information --region us-west-2 \
  --filters "Key=InstanceIds,Values=<instance-id>"
```

### AWS CLI Not Found

If you see errors about `aws` command not found, install the AWS CLI:

```bash
# macOS
brew install awscli

# Or download from: https://aws.amazon.com/cli/
```

### GNU Parallel Not Found

If you see errors about `parallel` or `sem` commands (used for concurrent S3 operations with rdev), install GNU parallel:

```bash
# macOS
brew install parallel

# Other systems: https://www.gnu.org/software/parallel/
```

### Tunnel Cleanup

Tunnels should clean up automatically via trap handlers ([mirror_env_data.sh#L42](mirror_env_data.sh#L42)), but if they don't:

```bash
# Find tunnel processes
ps aux | grep 'aws ssm start-session'

# Kill by PID
kill <PID>

# Or kill all SSM sessions
pkill -f 'aws ssm start-session'
```

### Check Tunnel Logs

If a tunnel fails to start, check the logs (log files created in [manage_db_tunnel.sh#L45](manage_db_tunnel.sh#L45)):

```bash
cat /tmp/tunnel-prod-5432.log
cat /tmp/tunnel-staging-5433.log
# Or for specific environment
cat /tmp/tunnel-<env>-<port>.log
```

## Advanced Usage

### Custom Port Numbers

If you need to run multiple mirror operations simultaneously (not recommended), use different ports:

```bash
# First mirror
make mirror_env_data SRC_ENV=prod DEST_ENV=staging SRC_PORT=5432 DEST_PORT=5433

# Second mirror (in another terminal)
make mirror_env_data SRC_ENV=prod DEST_ENV=dev SRC_PORT=5434 DEST_PORT=5435
```

### Manual Tunnel Management (Debug Mode)

If you need to debug tunnel issues, you can use the helper script directly ([manage_db_tunnel.sh](manage_db_tunnel.sh)):

```bash
# Start a tunnel
backend/scripts/manage_db_tunnel.sh start prod 5432

# Get the PID (printed to stdout)
# Example: 12345

# Stop the tunnel
backend/scripts/manage_db_tunnel.sh stop 12345
```

## Performance Notes

- **Database operations**: Dump/restore time depends on database size (typically 5-30 minutes)
- **S3 sync**: Time depends on data size (can take hours for large datasets)
  - Uses 16 parallel jobs for concurrent transfers ([mirror_s3_data.sh#L131](mirror_s3_data.sh#L131))
  - For rdev with COLLECTIONS, uses unlimited parallel jobs based on CPU cores ([mirror_s3_data.sh#L95](mirror_s3_data.sh#L95))
  - Uses AWS CLI with `--copy-props metadata-directive` to copy objects efficiently ([mirror_s3_data.sh#L25](mirror_s3_data.sh#L25))
  - Progress display is suppressed with `--no-progress` flag for cleaner output
- **Dry-run mode**: Much faster since database operations are skipped and S3 operations only check what would be changed without transferring data
- The operation may appear to hang during S3 sync - this is normal, be patient!

## Getting Help

For issues or questions:

1. Check the troubleshooting section above
2. Review tunnel logs in `/tmp/tunnel-*.log`
3. Check the Makefile comments for parameter details
4. Contact the platform team
