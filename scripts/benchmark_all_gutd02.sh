#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-8"

source "$SCRIPTDIR/bench_common.sh"

# other variables
export DATADIR=/data/SimulatedDatasets/Gut
export SIGS="$DATADIR/signal/PromethION_R10.4.1-seq2squiggle/d0.2"
export RES=$OUTDIR/benchmarks/gut_d2_$SIGNAL_LENGTH
export MANIFEST=$DATADIR/d0.2_manifest.tsv
export BASECALL_ADDRESS="ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock"
export BASECALL_CONFIG="dna_r10.4.1_e8.2_400bps_fast@v5.2.0||"
export BATCH_SIZE=4096
export DEBUG_LOG=$LOGDIR/debug.log
export INDIR=$TMPDIR/gutd2
export NUM_THREADS=16
export DATASET_LABEL=gut_d2

f0() { check_basecall_socket; }

f1() {
    require_signal_length 3200
    echo "Benchmarking mappy_rs for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    fn_idx_in=$INDIR/mm/gutd0.2.mmi \
    n_threads=$NUM_THREADS"

    run_benchmark readfish.plugins.mappy_rs "$pl_args" mappy_rs "$DATASET_LABEL"
}

f2() {
    require_signal_length 3200
    echo "Benchmarking metagraph query for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    method=query \
    input=$INDIR/mg/gutd0.2.dbg \
    annotator=$INDIR/mg/gutd0.2.column.annodbg \
    threads=$NUM_THREADS \
    num_top_labels=1 \
    discovery_fraction=0.1"

    run_benchmark pymetagraph "$pl_args" metagraph "$DATASET_LABEL"
}

f3() {
    require_signal_length 3200
    echo "Not Benchmarking metagraph align for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."
    return 0

    pl_args="\
    method=align \
    input=$INDIR/mg/gutd0.2.dbg \
    threads=$NUM_THREADS \
    seed_length=19 \
    max_alternative_alignments=1 \
    max_num_nodes_per_seq_char=5 \
    min_exact_match=0.5 \
    xdrop=50 \
    connect_anchors=false \
    extend_chains=false"

    run_benchmark pymetagraph "$pl_args" metagraph_align "$DATASET_LABEL" --max-batches 3
}

f4() {
    require_signal_length 3200
    echo "Benchmarking spumoni for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    ref=$INDIR/sp/gutd0.2 \
    threads=$NUM_THREADS \
    PML=true \
    minimizer_alphabet=true"

    run_benchmark pyspumoni "$pl_args" spumoni "$DATASET_LABEL"
}

f5() {
    require_signal_length 3200
    echo "Benchmarking readbouncer for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    target_files=$INDIR/rb/Refs_d0.2_Comm_1_k17.ibf \
    kmer_size=17 \
    threads=$NUM_THREADS"

    run_benchmark pyreadbouncer "$pl_args" readbouncer "$DATASET_LABEL"
}

f6() {
    require_signal_length 3200
    echo "Benchmarking collinearity for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    input=$INDIR/cl/gutd0.2.cidx \
    n_threads=$NUM_THREADS"

    run_benchmark pycollinearity "$pl_args" collinearity "$DATASET_LABEL"
}

f7() {
    require_signal_length 3200
    echo "Benchmarking rawhash for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    idx=$INDIR/rh/gutd0.2.ind \
    threads=$NUM_THREADS \
    x=faster"

    run_benchmark pyrawhash "$pl_args" rawhash "$DATASET_LABEL" --no-basecall --max-batches 10
}

f8() {
    require_signal_length 3200
    echo "Benchmarking sigmoni for Gut (d=0.2) with signal length $SIGNAL_LENGTH .."

    pl_args="\
    ref_prefix=$INDIR/sg/refs/ref \
    spumoni_path=$CODEDIR/spumoni/build/ \
    threads=$NUM_THREADS \
    multi=true \
    complexity=true"

    run_benchmark sigmoni "$pl_args" sigmoni "$DATASET_LABEL" --no-basecall --max-batches 100
}
