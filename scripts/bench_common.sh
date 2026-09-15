#!/bin/bash
#
# Shared helpers for benchmark_all_*.sh task files. Sourced via
# "$SCRIPTDIR/bench_common.sh" -- SCRIPTDIR is exported by task_runner.sh
# before it sources the task file, so this only works when invoked through
# task_runner.sh.
#

# require_signal_length [default]
# Validates SIGNAL_LENGTH is set and a positive multiple of 1600.
# - If valid: returns 0, SIGNAL_LENGTH untouched.
# - If invalid/unset and a default is given: sets SIGNAL_LENGTH to it, returns 0.
# - If invalid/unset and no default is given: returns 1 (caller should skip
#   the task, e.g. `require_signal_length || return 1`).
require_signal_length() {
    if [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); then
        return 0
    fi
    if [[ $# -ge 1 ]]; then
        SIGNAL_LENGTH=$1
        return 0
    fi
    return 1
}

# check_basecall_socket
# Manual connectivity check for the dorado basecall server socket (task f0 in
# every benchmark_all_*.sh).
check_basecall_socket() {
    # Does the socket exist and is it accessible from the calling script's environment?
    ls -la /var/lib/minknow/data/.dorado/dorado-basecall-server.sock

    # What user is the calling script running as?
    echo "Running as: $(whoami)"

    # Can you connect to the socket at all?
    python3 -c "
from pybasecall_client_lib.pyclient import PyBasecallClient
c = PyBasecallClient(
    address='ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock',
    config='dna_r10.4.1_e8.2_400bps_fast@v5.2.0||',
    priority=PyBasecallClient.high_priority,
    client_name='test',
)
print('connecting...')
c.connect()
print('connected')
c.disconnect()
"
}

# run_benchmark <plugin> <pl_args> <output_subdir> <dataset_label> [--no-basecall] [--max-batches N]
# The benchmark_aligner.py invocation shared by every benchmark task. Requires
# SIGS, RES, MANIFEST, BASECALL_ADDRESS, BASECALL_CONFIG, SIGNAL_LENGTH,
# BATCH_SIZE to already be set by the caller/task file.
run_benchmark() {
    local plugin="$1" pl_args="$2" outsubdir="$3" dataset_label="$4"; shift 4

    local basecall_args=(--basecall-address "$BASECALL_ADDRESS" --basecall-config "$BASECALL_CONFIG")
    local extra_args=()

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --no-basecall) basecall_args=(); shift ;;
            --max-batches) extra_args+=(--max-batches-per-file "$2"); shift 2 ;;
            *) echo "run_benchmark: unknown option '$1'" >&2; return 1 ;;
        esac
    done

    benchmark_aligner.py \
        --input "$SIGS/*.blow5" \
        --plugin "$plugin" \
        --plugin-args "$pl_args" \
        --output "$RES/$outsubdir/" \
        --manifest "$MANIFEST" \
        "${basecall_args[@]}" \
        --truncate-signals "$SIGNAL_LENGTH" \
        --batch-size "$BATCH_SIZE" \
        "${extra_args[@]}" \
        --dataset "$dataset_label"
}
