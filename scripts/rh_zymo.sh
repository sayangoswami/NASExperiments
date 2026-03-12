#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="0-1"

# other variables
export DATADIR=$EXPDIR/data/simulated/Zymo
export RAWHASH_DIR=~/Repos/RawHash/
export RAWHASH=$RAWHASH_DIR/rawhash2

f0() {
    measure $RAWHASH -t 32 -x viral -d $TMPDIR/zymo.ind -p \
        $RAWHASH_DIR/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model \
        $DATADIR/Refs1.fasta
}

f1() {
    measure $RAWHASH -x viral $TMPDIR/zymo.ind $DATADIR/signals/Sigs0_180.blow5 > $OUTDIR/zymo-rh.paf
    measure $RAWHASH -x viral $TMPDIR/zymo.ind $DATADIR/signals/Sigs1_180.blow5 >> $OUTDIR/zymo-rh.paf
}