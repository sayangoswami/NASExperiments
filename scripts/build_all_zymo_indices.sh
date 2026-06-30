#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-3"

# other variables
export DATADIR=/data/SimulatedDatasets/Zymo

f1() {
    echo "Building Collinearity index for Zymo dataset.."
    measure collinearity \
    --ref $DATADIR/Refs1.fasta \
    --idx $TMPDIR/Zymo
}

f2() {
    echo "Building Spumoni index for Zymo dataset.."
    spumoni build -r $DATADIR/Refs1.fasta -M -P -m -o $TMPDIR/zymo
}

f3() {
    echo "Building Metagraph index for Zymo dataset.."
    #  metagraph build and annotate
    metagraph build -v -p 16 -k 31 -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta
    metagraph annotate -v -p 16 -i $TMPDIR/mg_zymo.dbg --anno-header -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta 
}

