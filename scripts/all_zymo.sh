#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="0-2"

# other variables
export DATADIR=$EXPDIR/data/simulated/Zymo

f0() {
    spumoni build -r $DATADIR/Refs1.fasta -M -P -m -o $TMPDIR/zymo
}

f1() {
    spumoni run -r $TMPDIR/zymo -p $DATADIR/reads/Reads01_180.fasta -m -P -c -t 16
}

f2() {
    #  metagraph build and annotate
    metagraph build -v -p 16 -k 31 -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta
    metagraph annotate -v -p 16 -i $TMPDIR/mg_zymo.dbg --anno-header -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta 
}

f3() {
    # metagraph align
    metagraph align -p 16 \
        -i $TMPDIR/mg_zymo.dbg $DATADIR/reads/Reads01_180.fasta > $OUTDIR/zymo-mg.tsv
}