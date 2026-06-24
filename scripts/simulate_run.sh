#!/usr/bin/env bash
# set -x
set -euo pipefail

# --- Some globals ---
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXPDIR=`dirname $SCRIPTDIR`
CODEDIR=$EXPDIR/code
LOGDIR=$EXPDIR/logs
OUTDIR=$EXPDIR/results
TMPDIR=$EXPDIR/tmp

DATADIR=/data/SimulatedDatasets/Zymo/PromethION_R10.4.1-seq2squiggle/signals/

# --- Update INPUT_SIGNAL and CONFIG_TOML ---
# the following lines can be files or directories
INPUT_SIGNAL=$DATADIR
PROFILE="dna-r10-prom"
# the following is the path to Readfish's config.
CONFIG_TOML=/scratch/NASExperiments/configs/rf_rh_zymo.toml

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
MEM_LOG="$LOGDIR/readfish_mem_$TIMESTAMP.tsv"   # RSS log consumed by analyse_runs.py
SERVER_PID_FILE="$TMPDIR/server.pid"
MEM_POLLER_PID_FILE="$TMPDIR/mem_poller.pid"

# Memory polling interval in seconds. 1s is fine for a run that lasts minutes;
# increase to 5 if you are worried about overhead on a loaded machine.
MEM_POLL_INTERVAL=1

args=()
for f in "${INPUT_SIGNAL[@]}"; do
    args+=(--input "$f")
done

CLIENT_CMD="readfish targets --wait-for-ready 5 \
    --toml $CONFIG_TOML --port 50051 --device MN12345 \
    --log-file $CLIENT_LOG --experiment-name 'Readfish_$TIMESTAMP'"


# ---------------------------------------------------------------------------
# Memory poller
#
# Measures the RSS of the full Readfish process tree (main process + all
# children, recursively). This matters because Readfish's C++ .so plugins
# are loaded in child processes whose memory would be invisible if we only
# read /proc/<pid>/status.
#
# Per-process RSS is read from smaps_rollup (kernel ≥ 4.14, accurate) with
# a fallback to summing /proc/<pid>/smaps (older kernels, slower but correct).
# Using smaps rather than status/VmRSS means shared library pages are
# counted per-process rather than being collapsed, which slightly overstates
# total RSS when children share mappings — but for our purposes (comparing
# tools against each other under identical conditions) this is consistent and
# conservative, and it's the same methodology across all runs.
#
# Writes tab-separated (timestamp_s, rss_mb) rows to MEM_LOG, matching the
# format expected by analyse_runs.py's load_mem_log().
# ---------------------------------------------------------------------------

# Return all PIDs in the process tree rooted at $1 (including $1 itself).
get_process_tree() {
    local root=$1
    local pids=("$root")
    local children
    # pgrep -P lists direct children; we recurse until no more are found.
    children=$(pgrep -P "$root" 2>/dev/null || true)
    for child in $children; do
        pids+=( $(get_process_tree "$child") )
    done
    echo "${pids[@]}"
}

# Return total RSS in kB for a single PID, using smaps_rollup if available,
# falling back to summing the Rss: lines in /proc/<pid>/smaps.
rss_kb_for_pid() {
    local pid=$1
    if [[ -r /proc/$pid/smaps_rollup ]]; then
        awk '/^Rss:/{sum+=$2} END{print sum+0}' /proc/"$pid"/smaps_rollup 2>/dev/null
    elif [[ -r /proc/$pid/smaps ]]; then
        awk '/^Rss:/{sum+=$2} END{print sum+0}' /proc/"$pid"/smaps 2>/dev/null
    else
        echo 0
    fi
}

start_mem_poller() {
    local root_pid=$1
    echo -e "timestamp_s\trss_mb" > "$MEM_LOG"
    (
        while kill -0 "$root_pid" 2>/dev/null; do
            # Collect the full process tree at each sample point, since
            # child processes may be spawned or exit during the run.
            total_rss_kb=0
            for pid in $(get_process_tree "$root_pid"); do
                pid_rss=$(rss_kb_for_pid "$pid")
                total_rss_kb=$(( total_rss_kb + pid_rss ))
            done
            rss_mb=$(echo "scale=2; $total_rss_kb / 1024" | bc)
            echo -e "$(date +%s.%N)\t$rss_mb" >> "$MEM_LOG"
            sleep "$MEM_POLL_INTERVAL"
        done
    ) &
    MEM_POLLER_PID=$!
    echo "$MEM_POLLER_PID" > "$MEM_POLLER_PID_FILE"
    echo "Memory poller started (PID $MEM_POLLER_PID, tracking process tree of $root_pid every ${MEM_POLL_INTERVAL}s → $MEM_LOG)"
}

stop_mem_poller() {
    if [[ -f "$MEM_POLLER_PID_FILE" ]]; then
        MEM_POLLER_PID=$(<"$MEM_POLLER_PID_FILE")
        if kill -0 "$MEM_POLLER_PID" 2>/dev/null; then
            kill "$MEM_POLLER_PID"
            wait "$MEM_POLLER_PID" 2>/dev/null || true
        fi
        rm -f "$MEM_POLLER_PID_FILE"
    fi
}


# ---------------------------------------------------------------------------
# Cleanup handler for Ctrl+C or errors
# ---------------------------------------------------------------------------
cleanup() {
    echo "🧹 Cleaning up..."
    stop_mem_poller
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


# ---------------------------------------------------------------------------
# Dorado
# ---------------------------------------------------------------------------
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
    echo "Do you wish to continue without dorado? (y/n)"
    read -r answer
    if [[ "$answer" != "y" ]]; then
        exit 1
    fi
fi


# ---------------------------------------------------------------------------
# Start simulator (background)
# ---------------------------------------------------------------------------
echo "🚀 Starting server..."
mksimserver --profile "$PROFILE" --certs "$CERTS_DIR" "${args[@]}" >"$SERVER_LOG" 2>&1 &
SERVER_PID=$!
echo "$SERVER_PID" >"$SERVER_PID_FILE"
echo "Server PID: $SERVER_PID (logging to $SERVER_LOG)"

sleep 2


# ---------------------------------------------------------------------------
# Start Readfish (foreground) with memory polling
# We launch Readfish in the background just long enough to capture its PID,
# start the poller, then wait on it — so it still behaves as the foreground
# process from the user's perspective.
# ---------------------------------------------------------------------------
echo "💻 Starting client..."
$CLIENT_CMD &
CLIENT_PID=$!

# Start polling Readfish's RSS now that we have its PID
start_mem_poller "$CLIENT_PID"

# Wait for Readfish to finish — this blocks until it exits
wait "$CLIENT_PID"
CLIENT_EXIT=$?

# Stop the poller cleanly as soon as Readfish exits
stop_mem_poller

echo "Client finished (exit code $CLIENT_EXIT). Logs in $CLIENT_LOG"
echo "Memory log: $MEM_LOG"

# Print a quick peak-RSS summary directly to stdout
if [[ -f "$MEM_LOG" ]]; then
    peak=$(awk 'NR>1 {if($2>max) max=$2} END {printf "%.1f", max}' "$MEM_LOG")
    echo "Peak Readfish RSS: ${peak} MB"
fi


# ---------------------------------------------------------------------------
# Wait for simulator to finish
# ---------------------------------------------------------------------------
wait "$SERVER_PID" || true
echo "Server exited."