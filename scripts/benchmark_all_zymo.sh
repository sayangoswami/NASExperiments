#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-8"

source "$SCRIPTDIR/bench_common.sh"

# other variables
export DATADIR=/data/SimulatedDatasets/Zymo
export SIGS=$DATADIR/PromethION_R10.4.1-seq2squiggle/signals/
export RES=$OUTDIR/benchmarks/zymo_$SIGNAL_LENGTH
export MANIFEST=$DATADIR/manifest.tsv
export BASECALL_ADDRESS="ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock"
export BASECALL_CONFIG="dna_r10.4.1_e8.2_400bps_fast@v5.2.0||"
export BATCH_SIZE=4096
export DEBUG_LOG=$LOGDIR/debug.log
export INDIR=$TMPDIR/zymo
export NUM_THREADS=16
export DATASET_LABEL=zymo

f0() { check_basecall_socket; }

f1() {
    require_signal_length || return 1
    echo "Benchmarking mappy_rs for Zymo.."

    pl_args="\
    fn_idx_in=$INDIR/mm/zymo.mmi \
    n_threads=$NUM_THREADS"

    run_benchmark readfish.plugins.mappy_rs "$pl_args" mappy_rs "$DATASET_LABEL"
}

f2() {
    require_signal_length || return 1
    echo "Benchmarking metagraph query for Zymo.."

    pl_args="\
    method=query \
    input=$INDIR/mg/zymo.dbg \
    annotator=$INDIR/mg/zymo.column.annodbg \
    threads=$NUM_THREADS \
    num_top_labels=1 \
    discovery_fraction=0.1"

    run_benchmark pymetagraph "$pl_args" metagraph "$DATASET_LABEL"
}

f3() {
    require_signal_length || return 1
    echo "Not Benchmarking metagraph align for Zymo.."
    return 0

    pl_args="\
    method=align \
    input=$INDIR/mg/zymo.dbg \
    threads=$NUM_THREADS \
    seed_length=19 \
    max_alternative_alignments=1 \
    max_num_nodes_per_seq_char=5 \
    min_exact_match=0.5 \
    connect_anchors=false \
    extend_chains=false"

    run_benchmark pymetagraph "$pl_args" metagraph_align "$DATASET_LABEL"
}

f4() {
    require_signal_length || return 1
    echo "Benchmarking spumoni for Zymo.."

    pl_args="\
    ref=$INDIR/sp/zymo \
    threads=$NUM_THREADS \
    PML=true \
    minimizer_alphabet=true"

    run_benchmark pyspumoni "$pl_args" spumoni "$DATASET_LABEL"
}

f5() {
    require_signal_length || return 1
    echo "Benchmarking readbouncer for Zymo.."

    pl_args="\
    kmer_size=17 \
    target_files=$INDIR/rb/Refs1.ibf \
    exp_seq_error_rate=0.05 \
    threads=$NUM_THREADS"

    run_benchmark pyreadbouncer "$pl_args" readbouncer "$DATASET_LABEL"
}

f6() {
    require_signal_length || return 1
    echo "Benchmarking collinearity for Zymo.."

    pl_args="\
    input=$INDIR/cl/zymo.cidx \
    n_threads=$NUM_THREADS"

    run_benchmark pycollinearity "$pl_args" collinearity "$DATASET_LABEL"
}

f7() {
    require_signal_length || return 1
    echo "Benchmarking rawhash for Zymo.."

    pl_args="\
    idx=$INDIR/rh/zymo.ind \
    threads=$NUM_THREADS \
    x=viral"

    run_benchmark pyrawhash "$pl_args" rawhash "$DATASET_LABEL" --no-basecall
}

f8() {
    require_signal_length || return 1
    echo "Benchmarking sigmoni for Zymo.."

    pl_args="\
    ref_prefix=$INDIR/sg/refs/ref \
    spumoni_path=$CODEDIR/spumoni/build/ \
    threads=$NUM_THREADS"

    run_benchmark sigmoni "$pl_args" sigmoni "$DATASET_LABEL" --no-basecall
}
