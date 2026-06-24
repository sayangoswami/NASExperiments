#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="0-1"

# other variables
export DATADIR=/data/SimulatedDatasets/Zymo
export REF=$DATADIR/Refs1.fasta
export SIGS=$DATADIR/PromethION_R10.4.1-seq2squiggle/signals/
export RAWHASH_DIR=$CODEDIR/rawhash2

f0() {
    measure rawhash2 --r10 -t 32 -d $TMPDIR/zymo.ind -p \
        $RAWHASH_DIR/extern/local_kmer_models/uncalled_r1041_model_only_means.txt \
        $DATADIR/Refs1.fasta
}

f1() {
    benchmark_aligner.py \
        --input $SIGS/ \
        --plugin pyrawhash \
        --plugin-args "idx=$TMPDIR/zymo.ind threads=16 x=viral" \
        --output $RESDIR/rawhash_zymo/
}

f2() {
    measure rawhash2 -x viral $TMPDIR/zymo.ind $SIGS/Sigs0_1000.blow5 > $OUTDIR/zymo-rh.paf
    measure rawhash2 -x viral $TMPDIR/zymo.ind $SIGS/Sigs1_1000.blow5 >> $OUTDIR/zymo-rh.paf
}