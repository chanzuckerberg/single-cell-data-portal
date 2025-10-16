#!/usr/bin/env bash

# Helper script to manage SSH tunnels for database access
# Usage: 
#   ./manage_db_tunnel.sh start <env> <port>  - Start tunnel and return PID
#   ./manage_db_tunnel.sh stop <pid>           - Stop tunnel by PID

set -e

SCRIPTS_DIR=$(dirname "$0")
source "${SCRIPTS_DIR}/_mirror_utils.sh"

ACTION=$1
DEPLOYMENT_STAGE=$2
PORT=$3

if [[ "$ACTION" == "start" ]]; then
    # Set AWS profile and get instance environment
    set_aws_profile "$DEPLOYMENT_STAGE"
    INSTANCE_ENV=$(get_instance_env "$DEPLOYMENT_STAGE")

    # Determine cluster name
    if [[ $DEPLOYMENT_STAGE == 'rdev' ]]; then
        CLUSTER_NAME="dataportal-rdev-happy"
    else
        CLUSTER_NAME="corpora-${DEPLOYMENT_STAGE}-corpora-api"
    fi

    # Get RDS endpoint
    endpoint=$(aws rds describe-db-cluster-endpoints --region ${AWS_REGION} --db-cluster-identifier ${CLUSTER_NAME} | jq -r '.DBClusterEndpoints[] | select(.EndpointType | contains("WRITER")) | .Endpoint')
    
    # Get EC2 instance ID
    instance_id=$(aws ec2 describe-instances --region ${AWS_REGION} --filters "Name=tag:Name,Values=dp-${INSTANCE_ENV}-happy" --query "Reservations[*].Instances[*].InstanceId" --output text)

    echo "Starting tunnel for $DEPLOYMENT_STAGE on port $PORT..." >&2
    echo "  Endpoint: $endpoint" >&2
    echo "  Instance: $instance_id" >&2

    # Start SSH tunnel in background and capture PID
    nohup aws ssm start-session \
        --region ${AWS_REGION} \
        --target ${instance_id} \
        --document-name AWS-StartPortForwardingSessionToRemoteHost \
        --parameters "{\"portNumber\":[\"5432\"],\"localPortNumber\":[\"${PORT}\"],\"host\":[\"${endpoint}\"]}" \
        > "/tmp/tunnel-${DEPLOYMENT_STAGE}-${PORT}.log" 2>&1 &
    
    tunnel_pid=$!
    
    # Wait a few seconds for tunnel to establish
    echo "Waiting for tunnel to establish..." >&2
    sleep 5
    
    # Check if process is still running
    if ps -p $tunnel_pid > /dev/null 2>&1; then
        echo "Tunnel started successfully with PID: $tunnel_pid" >&2
        echo $tunnel_pid
    else
        echo "ERROR: Tunnel failed to start. Check /tmp/tunnel-${DEPLOYMENT_STAGE}-${PORT}.log" >&2
        exit 1
    fi

elif [[ "$ACTION" == "stop" ]]; then
    PID=$2
    if [[ -z "$PID" ]]; then
        echo "ERROR: No PID provided to stop" >&2
        exit 1
    fi
    
    echo "Stopping tunnel with PID: $PID" >&2
    
    # Kill the process and its children
    if ps -p $PID > /dev/null 2>&1; then
        # Get all child processes
        pkill -P $PID 2>/dev/null || true
        sleep 1
        kill $PID 2>/dev/null || true
        echo "Tunnel stopped" >&2
    else
        echo "Process $PID not found (may have already terminated)" >&2
    fi

else
    echo "Usage: $0 {start|stop} <args...>" >&2
    echo "  start <env> <port>  - Start tunnel for environment on specified port" >&2
    echo "  stop <pid>          - Stop tunnel by PID" >&2
    exit 1
fi

