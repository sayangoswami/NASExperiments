#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-7"

# other variables
export DATADIR=/data/SimulatedDatasets/Gut/
export RAWHASH_DIR=$CODEDIR/rawhash2
export OUTDIR=$TMPDIR/gutd2
export REF=$DATADIR/Refs_d0.2_Comm_1.fa

f1() {
    mkdir -p $OUTDIR/cl
    echo "Building Collinearity index for Gut (d=0.2).."
    measure Collinearity \
    --ref $REF \
    --idx $OUTDIR/cl/gutd0.2 \
    --bw 1024
}

f2() {
    mkdir -p $OUTDIR/sp
    echo "Building Spumoni index for Gut (d=0.2).."
    measure spumoni build -r $REF -M -P -m -o $OUTDIR/sp/gutd0.2
}

f3() {
    mkdir -p $OUTDIR/mg
    echo "Building Metagraph index for Gut (d=0.2).."
    #  metagraph build and annotate
    measure metagraph build -v -p 16 -k 21 -o $OUTDIR/mg/gutd0.2 $REF
    measure metagraph annotate -v -p 16 -i $OUTDIR/mg/gutd0.2.dbg --anno-header -o $OUTDIR/mg/gutd0.2 $REF 
}

f4() {
    mkdir -p $OUTDIR/mm
    echo "Building Minimap2 index for Gut (d=0.2).."
    measure minimap2 -x map-ont -d $OUTDIR/mm/gutd0.2.mmi $REF
}

f5() {
    mkdir -p $OUTDIR/rb
    echo "Building readBouncer index for Gut (d=0.2).."
    local config_file
    config_file="$(mktemp --suffix=.toml)"

    # cleanup fires when this function returns, however it returns
    trap 'rm -f "$config_file"' RETURN

    cat > "$config_file" <<EOF
usage               = "build"
output_directory    = "$OUTDIR/rb"
log_directory       = "$LOGDIR"

[IBF]
threads       = 8
target_files  = ["$REF"]
EOF

    measure ReadBouncer --config "$config_file"
}

f6() {
    mkdir -p $OUTDIR/rh
    echo "Building RawHash index for Gut (d=0.2).."
    measure rawhash2 --r10 -t 32 -d $OUTDIR/rh/gutd0.2.ind -p \
        $RAWHASH_DIR/extern/local_kmer_models/uncalled_r1041_model_only_means.txt \
        $REF
}

f7() {
    mkdir -p $OUTDIR/sg
    echo "Building Sigmoni index for Gut (d=0.2).."
    measure sigmoni-index -p $REF \
        -o $OUTDIR/sg --spumoni-path $SPUMONI_BUILD_DIR \
        --poremodel $CODEDIR/sigmoni/poremodel/model_r1041_400bps_dm_it2.tsv
}