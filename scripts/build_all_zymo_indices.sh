#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-7"

# other variables
export DATADIR=/data/SimulatedDatasets/Zymo
export RAWHASH_DIR=$CODEDIR/rawhash2

f1() {
    echo "Building Collinearity index for Zymo.."
    measure Collinearity \
    --ref $DATADIR/Refs1.fasta \
    --idx $TMPDIR/Zymo
}

f2() {
    echo "Building Spumoni index for Zymo.."
    measure spumoni build -r $DATADIR/Refs1.fasta -M -P -m -o $TMPDIR/zymo
}

f3() {
    echo "Building Metagraph index for Zymo.."
    #  metagraph build and annotate
    measure metagraph build -v -p 16 -k 31 -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta
    measure metagraph annotate -v -p 16 -i $TMPDIR/mg_zymo.dbg --anno-header -o $TMPDIR/mg_zymo $DATADIR/Refs1.fasta 
}

f4() {
    echo "Building Minimap2 index for Zymo.."
    measure minimap2 -x map-ont -d $TMPDIR/zymo.mmi $DATADIR/Refs1.fasta
}

f5() {
    echo "Building readBouncer index for Zymo.."
    local config_file
    config_file="$(mktemp --suffix=.toml)"

    # cleanup fires when this function returns, however it returns
    trap 'rm -f "$config_file"' RETURN

    cat > "$config_file" <<EOF
usage               = "build"
output_directory    = "$TMPDIR"
log_directory       = "$LOGDIR"

[IBF]
threads       = 8
target_files  = ["$DATADIR/Refs1.fasta"]
EOF

    measure ReadBouncer --config "$config_file"
}

f6() {
    echo "Building RawHash index for Zymo.."
    measure rawhash2 --r10 -t 32 -d $TMPDIR/zymo.ind -p \
        $RAWHASH_DIR/extern/local_kmer_models/uncalled_r1041_model_only_means.txt \
        $DATADIR/Refs1.fasta
}

f7() {
    echo "Building Sigmoni index for Zymo.."
    measure sigmoni-index -p $DATADIR/Refs1.fasta \
        -o $TMPDIR --spumoni-path $SPUMONI_BUILD_DIR
}