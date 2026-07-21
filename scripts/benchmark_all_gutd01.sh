#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-7"

# other variables
export SIGNAL_LENGTH=3200
export DATADIR=/data/SimulatedDatasets/Gut
export SIGS="$DATADIR/signal/PromethION_R10.4.1-seq2squiggle/d0.1"
export RES=$RESDIR/benchmarks/gut_d1_$SIGNAL_LENGTH
export MANIFEST=$DATADIR/d0.1_manifest.tsv
export BASECALL_ADDRESS="ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock"
export BASECALL_CONFIG="dna_r10.4.1_e8.2_400bps_fast@v5.2.0||"
export BATCH_SIZE=4096
export DEBUG_LOG=$LOGDIR/debug.log
export INDIR=$TMPDIR/gutd1

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
    echo "Benchmarking mappy_rs for Gut (d=0.1).."
    
    pl_args="\
    fn_idx_in=$INDIR/mm/gutd0.1.mmi \
    n_threads=16"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin readfish.plugins.mappy_rs \
    --plugin-args "$pl_args" \
    --output $RES/mappy_rs/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f2() {
    echo "Benchmarking metagraph query for Gut (d=0.1).."
    
    pl_args="\
    method=query \
    input=$INDIR/mg/gutd0.1.dbg \
    annotator=$INDIR/mg/gutd0.1.column.annodbg \
    threads=16 \
    num_top_labels=1 \
    discovery_fraction=0.1"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pymetagraph \
    --plugin-args "$pl_args" \
    --output $RES/metagraph_query/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f3() {
    echo "Benchmarking metagraph align for Gut (d=0.1).."

    pl_args="\
    method=align \
    input=$INDIR/mg/gutd0.1.dbg \
    annotator=$INDIR/mg/gutd0.1.column.annodbg \
    threads=16 \
    seed_length=21 \
    max_alternative_alignments=1 \
    max_num_nodes_per_seq_char=10 \
    min_exact_match=0.4 \
    connect_anchors=true \
    extend_chains=true"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pymetagraph \
    --plugin-args "$pl_args" \
    --output $RES/metagraph_align/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f4() {
    echo "Benchmarking spumoni for Gut (d=0.1).."
    
    pl_args="\
    ref=$INDIR/sp/gutd0.1 \
    threads=16 \
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
    --batch-size $BATCH_SIZE
}

f5() {
    echo "Benchmarking readbouncer for Gut (d=0.1).."
    
    pl_args="\
    target_files=$INDIR/rb/Refs_d0.1_Comm_1.ibf \
    threads=16"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pyreadbouncer \
    --plugin-args "$pl_args" \
    --output $RES/readbouncer/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f6() {
    echo "Benchmarking collinearity for Gut (d=0.1).."
    
    pl_args="\
    input=$INDIR/cl/gutd0.1.cidx \
    bw=1024 \
    n_threads=16"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pycollinearity \
    --plugin-args "$pl_args" \
    --output $RES/collinearity/ \
    --manifest $MANIFEST \
    --basecall-address $BASECALL_ADDRESS \
    --basecall-config  $BASECALL_CONFIG \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f7() {
    echo "Benchmarking rawhash for Gut (d=0.1).."
    
    pl_args="\
    idx=$INDIR/rh/gutd0.1.ind \
    threads=16 \
    x=faster"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin pyrawhash \
    --plugin-args "$pl_args" \
    --output $RES/rawhash/ \
    --manifest $MANIFEST \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}

f8() {
    echo "Benchmarking sigmoni for Gut (d=0.1).."
    
    pl_args="\
    ref_prefix=$INDIR/sg/refs/ref \
    spumoni_path=$CODEDIR/spumoni/build/ \
    threads=16"
    
    benchmark_aligner.py \
    --input "$SIGS/*.blow5" \
    --plugin sigmoni \
    --plugin-args "$pl_args" \
    --output $RES/sigmoni/ \
    --manifest $MANIFEST \
    --truncate-signals $SIGNAL_LENGTH \
    --batch-size $BATCH_SIZE
}