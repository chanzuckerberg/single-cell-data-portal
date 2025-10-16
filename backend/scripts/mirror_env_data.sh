#!/usr/bin/env bash

# Mirror S3 and RDS (postgres db) data from a source deployment environment to a destination environment.
# Automatically manages SSH tunnels for RDS operations.
#
# Usage: See 'make mirror_env_data' in the Makefile or MIRROR_ENV_DATA_USAGE.md for details.
#
# THIS IS DESTRUCTIVE for the destination env! The source env will never be modified,
# but the destination env's data will be completely replaced.

set -e

SCRIPTS_DIR=$(dirname "$0")
source "${SCRIPTS_DIR}/_mirror_utils.sh"

# Source the env setup to get SRC_ENV, DEST_ENV, etc.
source "${SCRIPTS_DIR}/set_src_dest_envs.sh"

if [[ -n "$DRY_RUN" ]]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                     🔍 DRY-RUN MODE 🔍                         ║"
    echo "║                                                                ║"
    echo "║  Testing connections and showing what would be changed.       ║"
    echo "║  No data will be modified.                                    ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
fi

log "Starting environment data mirror: $SRC_ENV → $DEST_ENV"

# ============================================================
# PHASE 1: Mirror RDS Data
# ============================================================

log "Phase 1: RDS Database Mirroring"

SRC_TUNNEL_PID=""
DEST_TUNNEL_PID=""

# Cleanup function to stop tunnels
cleanup_tunnels() {
    if [[ -n "$SRC_TUNNEL_PID" ]]; then
        log "Stopping source tunnel (PID: $SRC_TUNNEL_PID)"
        $SCRIPTS_DIR/manage_db_tunnel.sh stop $SRC_TUNNEL_PID || true
    fi
    
    if [[ -n "$DEST_TUNNEL_PID" ]]; then
        log "Stopping destination tunnel (PID: $DEST_TUNNEL_PID)"
        $SCRIPTS_DIR/manage_db_tunnel.sh stop $DEST_TUNNEL_PID || true
    fi
}

# Set up trap to cleanup tunnels on exit
trap cleanup_tunnels EXIT INT TERM

# Start SSH tunnels for RDS access
log "Starting SSH tunnels for database access..."
SRC_TUNNEL_PID=$($SCRIPTS_DIR/manage_db_tunnel.sh start $SRC_ENV $SRC_PORT)
DEST_TUNNEL_PID=$($SCRIPTS_DIR/manage_db_tunnel.sh start $DEST_ENV $DEST_PORT)

log "Tunnels ready: src=$SRC_PORT (PID $SRC_TUNNEL_PID), dest=$DEST_PORT (PID $DEST_TUNNEL_PID)"
sleep 3  # Let tunnels stabilize

# Run RDS mirroring if not in dry-run mode
export NO_PROMPT=1  # Don't prompt again
if [[ -z "$DRY_RUN" ]]; then
    $SCRIPTS_DIR/mirror_rds_data.sh $@
fi

# ============================================================
# PHASE 2: Mirror S3 Data
# ============================================================

log "Phase 2: S3 Data Mirroring"

if [[ $DEST_ENV != 'rdev' ]]; then
  $SCRIPTS_DIR/mirror_s3_data.sh $@
elif [[ -n "$COLLECTIONS" ]]; then  # For rdev
  collections_array=($(tr ',' '\n' <<< $COLLECTIONS))
  id_regex='[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}'
  for uuid in "${collections_array[@]}"; do
    if [[ -z $(grep -E $id_regex <<< $uuid) ]]; then
      error_exit "$uuid is not a valid uuid"
    fi
  done
  $SCRIPTS_DIR/mirror_s3_data.sh $@
else
  log "DEST_ENV is set to rdev -- will NOT copy s3 assets"
fi

log "✅ Mirror operation completed successfully!"

# Clean up tunnels
log "mirroring complete. Closing SSH tunnels..."
cleanup_tunnels
trap - EXIT INT TERM  # Remove trap since we're handling cleanup manually
SRC_TUNNEL_PID=""
DEST_TUNNEL_PID=""
