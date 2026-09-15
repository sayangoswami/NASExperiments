#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-8"

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

f0() {
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

f1() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking mappy_rs for Zymo.."
    
    pl_args="\
    fn_idx_in=$INDIR/mm/zymo.mmi \
    n_threads=$NUM_THREADS"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin readfish.plugins.mappy_rs \
    --plugin-args "$pl_args" \
    --output $RES/mappy_rs/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f2() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking metagraph query for Zymo.."
    
    pl_args="\
    method=query \
    input=$INDIR/mg/zymo.dbg \
    annotator=$INDIR/mg/zymo.column.annodbg \
    threads=$NUM_THREADS \
    num_top_labels=1 \
    discovery_fraction=0.1"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pymetagraph \
    --plugin-args "$pl_args" \
    --output $RES/metagraph/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f3() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
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
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pymetagraph \
    --plugin-args "$pl_args" \
    --output $RES/metagraph_align/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f4() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking spumoni for Zymo.."
    
    pl_args="\
    ref=$INDIR/sp/zymo \
    threads=$NUM_THREADS \
    PML=true \
    minimizer_alphabet=true"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pyspumoni \
    --plugin-args "$pl_args" \
    --output $RES/spumoni/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f5() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking readbouncer for Zymo.."
    
    pl_args="\
    kmer_size     = 17
    target_files=$INDIR/rb/Refs1.ibf \
    exp_seq_error_rate=0.05 \
    threads=$NUM_THREADS"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pyreadbouncer \
    --plugin-args "$pl_args" \
    --output $RES/readbouncer/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f6() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking collinearity for Zymo.."
    
    pl_args="\
    input=$INDIR/cl/zymo.cidx \
    n_threads=$NUM_THREADS"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pycollinearity \
    --plugin-args "$pl_args" \
    --output $RES/collinearity/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f7() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking rawhash for Zymo.."
    
    pl_args="\
    idx=$INDIR/rh/zymo.ind \
    threads=$NUM_THREADS \
    x=viral"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pyrawhash \
    --plugin-args "$pl_args" \
    --output $RES/rawhash/ \
    --manifest $MANIFEST \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}

f8() {
    { [[ -v SIGNAL_LENGTH ]] && (( SIGNAL_LENGTH >= 1600 && SIGNAL_LENGTH % 1600 == 0 )); } || return 1
    echo "Benchmarking sigmoni for Zymo.."
    
    pl_args="\
    ref_prefix=$INDIR/sg/refs/ref \
    spumoni_path=$CODEDIR/spumoni/build/ \
    threads=$NUM_THREADS"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin sigmoni \
    --plugin-args "$pl_args" \
    --output $RES/sigmoni/ \
    --manifest $MANIFEST \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE \
    --dataset zymo
}