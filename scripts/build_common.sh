#!/bin/bash
#
# Shared helpers for build_all_*.sh task files. Sourced via
# "$SCRIPTDIR/build_common.sh" -- SCRIPTDIR is exported by task_runner.sh
# before it sources the task file, so this only works when invoked through
# task_runner.sh.
#

# build_start <output_subdir> <message>
# Common preamble for every build task: waits briefly (letting a previous
# task's background work settle), creates the tool's output subdirectory
# under $OUTDIR, and logs the given message.
build_start() {
    local subdir="$1" message="$2"
    sleep 2
    mkdir -p "$OUTDIR/$subdir"
    echo "$message"
}

# build_readbouncer_index <kmer_size>
# Writes a temporary ReadBouncer build config (removed when the calling
# function returns) and runs the build. Requires OUTDIR, LOGDIR, NUM_THREADS,
# REF to already be set.
build_readbouncer_index() {
    local kmer_size="$1"
    local config_file
    config_file="$(mktemp --suffix=.toml)"

    # cleanup fires when this function returns, however it returns
    trap 'rm -f "$config_file"' RETURN

    cat > "$config_file" <<EOF
usage               = "build"
output_directory    = "$OUTDIR/rb"
log_directory       = "$LOGDIR"

[IBF]
kmer_size     = $kmer_size
threads       = $NUM_THREADS
target_files  = ["$REF"]
EOF

    measure ReadBouncer --config "$config_file"
}
