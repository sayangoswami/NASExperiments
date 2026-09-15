#!/bin/bash

# --- List of all tasks to run by default ---
# This variable is used by the main script when no numbers are specified.
export ALL_TASKS="1-7"

# other variables
export DATADIR=/data/SimulatedDatasets/Zymo
export RAWHASH_DIR=$CODEDIR/rawhash2
export OUTDIR=$TMPDIR/zymo
export NUM_THREADS=32
export OMP_NUM_THREADS=$NUM_THREADS

f1() {
    sleep 2
    mkdir -p $OUTDIR/cl
    echo "Building Collinearity index for Zymo.."
    measure Collinearity \
    --ref $DATADIR/Refs1.fasta \
    --idx $OUTDIR/cl/zymo \
    --bw 1024 --sparse --n_threads $NUM_THREADS
}

f2() {
    sleep 2
    mkdir -p $OUTDIR/sp
    echo "Building Spumoni index for Zymo.."
    measure spumoni build -r $DATADIR/Refs1.fasta -M -P -m -o $OUTDIR/sp/zymo
}

f3() {
    sleep 2
    mkdir -p $OUTDIR/mg
    echo "Building Metagraph index for Zymo.."
    #  metagraph build and annotate
    measure metagraph build --mode canonical -v -p $NUM_THREADS -k 19 -o $OUTDIR/mg/zymo $DATADIR/Refs1.fasta
    measure metagraph annotate -v -p $NUM_THREADS -i $OUTDIR/mg/zymo.dbg --anno-header -o $OUTDIR/mg/zymo $DATADIR/Refs1.fasta 
}

f4() {
    sleep 2
    mkdir -p $OUTDIR/mm
    echo "Building Minimap2 index for Zymo.."
    measure minimap2 -t $NUM_THREADS -x map-ont -d $OUTDIR/mm/zymo.mmi $DATADIR/Refs1.fasta
}

f5() {
    sleep 2
    mkdir -p $OUTDIR/rb
    echo "Building readBouncer index for Zymo.."
    local config_file
    config_file="$(mktemp --suffix=.toml)"

    # cleanup fires when this function returns, however it returns
    trap 'rm -f "$config_file"' RETURN

    cat > "$config_file" <<EOF
usage               = "build"
output_directory    = "$OUTDIR/rb"
log_directory       = "$LOGDIR"

[IBF]
kmer_size     = 17
threads       = $NUM_THREADS
target_files  = ["$DATADIR/Refs1.fasta"]
EOF

    measure ReadBouncer --config "$config_file"
}

f6() {
    sleep 2
    mkdir -p $OUTDIR/rh
    echo "Building RawHash index for Zymo.."
    measure rawhash2 --r10 -t $NUM_THREADS -d $OUTDIR/rh/zymo.ind -p \
        $RAWHASH_DIR/extern/local_kmer_models/uncalled_r1041_model_only_means.txt \
        $DATADIR/Refs1.fasta
}

f7() {
    sleep 2
    mkdir -p $OUTDIR/sg
    echo "Building Sigmoni index for Zymo.."
    measure sigmoni-index -p $DATADIR/Refs1.fasta \
        -o $OUTDIR/sg --spumoni-path $SPUMONI_BUILD_DIR \
        --poremodel $CODEDIR/sigmoni/poremodel/model_r1041_400bps_dm_it2.tsv
}