#!/usr/bin/env bash

# Shared utilities for mirror scripts
# This file should be sourced by other mirror scripts

# AWS region used for all operations
export AWS_REGION="us-west-2"

# Get the appropriate AWS profile for a given environment
# Usage: get_aws_profile <env>
# Returns: The AWS profile name (e.g., "single-cell-prod", "single-cell-dev")
get_aws_profile() {
    local env=$1
    
    if [[ "$env" == "staging" || "$env" == "dev" || "$env" == "rdev" ]]; then
        echo "single-cell-dev"
    else
        echo "single-cell-${env}"
    fi
}

# Set AWS profile for a given environment
# Usage: set_aws_profile <env>
set_aws_profile() {
    local env=$1
    export AWS_PROFILE=$(get_aws_profile "$env")
}

# Get the instance environment name (differs from DEPLOYMENT_STAGE for staging)
# Usage: get_instance_env <env>
# Returns: Instance tag name (e.g., "stage" for staging, "dev" for dev)
get_instance_env() {
    local env=$1
    
    if [[ "$env" == "staging" ]]; then
        echo "stage"
    else
        echo "$env"
    fi
}

# Get database name and user for an environment
# Usage: get_db_credentials <env> <stack>
# Outputs: DB_NAME and DB_USER as space-separated values
get_db_credentials() {
    local env=$1
    local stack=$2
    
    if [[ "$env" == "rdev" ]]; then
        echo "/${stack} dataportal"
    else
        echo "corpora_${env} corpora_${env}"
    fi
}

# Log a message with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Log an error message and exit
error_exit() {
    echo "ERROR: $*" >&2
    exit 1
}

