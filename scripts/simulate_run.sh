#!/usr/bin/env bash
# set -x
set -euo pipefail

# --- Update INPUT_SIGNAL and CONFIG_TOML ---
# the following lines can be files or directories
INPUT_SIGNAL=(
    "data/simulated/Zymo/signals/Sigs0_450.blow5"
    "data/simulated/Zymo/signals/Sigs1_450.blow5"
)
# the folloing is the path to Readfish's config.
CONFIG_TOML=/scratch/NASExperiments/configs/rf_sp_zymo.toml

# --- Some globals ---
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXPDIR=`dirname $SCRIPTDIR`
CODEDIR=$EXPDIR/code
LOGDIR=$EXPDIR/logs
OUTDIR=$EXPDIR/results
TMPDIR=$EXPDIR/tmp

export PYTHONUNBUFFERED=1
export MINKNOW_API_USE_LOCAL_TOKEN="no" # to avoid printing a debug error message (non-fatal)
export MINKNOW_SIMULATOR="true"
CERTS_DIR=$CODEDIR/MinknoApiSimulator/certs
export MINKNOW_TRUSTED_CA=$CERTS_DIR/server.pem
export MINKNOW_API_CLIENT_CERTIFICATE_CHAIN=$CERTS_DIR/client.pem
export MINKNOW_API_CLIENT_KEY=$CERTS_DIR/client.key
export TIMESTAMP=$(date +"%b%d_%H%M%S")

# Configuration
SERVER_LOG="$LOGDIR/server_$TIMESTAMP.log"
CLIENT_LOG="$LOGDIR/readfish_$TIMESTAMP.log"
SERVER_PID_FILE="$TMPDIR/server.pid"

args=()
for f in "${INPUT_SIGNAL[@]}"; do
    args+=(--input "$f")
done

CLIENT_CMD="readfish targets --wait-for-ready 5 \
    --toml $CONFIG_TOML --port 50051 --device MN12345 \
    --log-file $CLIENT_LOG --experiment-name 'Readfish_$TIMESTAMP'"


# Cleanup handler for Ctrl+C or errors
cleanup() {
    echo "🧹 Cleaning up..."
    if [[ -f "$SERVER_PID_FILE" ]]; then
        SERVER_PID=$(<"$SERVER_PID_FILE")
        if kill -0 "$SERVER_PID" 2>/dev/null; then
            echo "Stopping server (PID $SERVER_PID)..."
            kill "$SERVER_PID"
            wait "$SERVER_PID" 2>/dev/null || true
        fi
        rm -f "$SERVER_PID_FILE"
    fi
}
trap cleanup EXIT INT TERM

# if the system dorado server is not running, and I have a local dorado, run it in the background
LOCAL_DORADO="$CODEDIR/ont-dorado-server/bin/dorado_basecall_server"

if pgrep -f -- "dorado_basecall_server" > /dev/null; then
    echo "System dorado server is already running"
elif [[ -x "$LOCAL_DORADO" ]]; then
    echo "Dorado is not running, starting local dorado server..."
    "$LOCAL_DORADO" --log_path $LOGDIR/dorado \
        --ipc_threads 3 --port /tmp/.guppy/5555 \
        --dorado_download_path $TMPDIR/dorado-models --device auto 2>&1 &
    DORADO_PID=$!
    sleep 2
    trap "kill $DORADO_PID; cleanup" EXIT INT TERM
else
    echo "System dorado server is not running, and local binary is missing: $LOCAL_DORADO"
    echo " Do you wish to continue without dorado? (y/n)"
    read -r answer
    if [[ "$answer" != "y" ]]; then
        exit 1
    fi
fi

# Start server in background
echo "🚀 Starting server..."
mksimserver --certs "$CERTS_DIR" "${args[@]}" >"$SERVER_LOG" 2>&1 &
SERVER_PID=$!
echo "$SERVER_PID" >"$SERVER_PID_FILE"
echo "Server PID: $SERVER_PID (logging to $SERVER_LOG)"

# Give the server a moment to start
sleep 2

# Start client (foreground)
echo "💻 Starting client..."
$CLIENT_CMD
echo "Client finished. Logs in $CLIENT_LOG"

# Wait for server to finish (if it exits)
wait "$SERVER_PID" || true
echo "Server exited."
