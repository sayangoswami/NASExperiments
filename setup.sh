#!/bin/bash

export EXPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $EXPDIR

# Create directory structure
mkdir -p code logs out results tmp

# setup minimap2
cd $EXPDIR/code
[ -d minimap2 ] || git clone https://github.com/lh3/minimap2.git
cd minimap2
make -j8

# setup collinearity
cd $EXPDIR/code
[ -d collinearity ] || git clone --recursive https://github.com/ratschlab/collinearity.git
cd collinearity

# for python bindings
pip install .

# to build from source
mkdir -p build && cd build
cmake ..
make -j 8

# build metagraph from sg/readfish-bindings branch
# From https://metagraph.ethz.ch/static/docs/installation.html#install-from-source
cd $EXPDIR/code
[ -d metagraph ] || git clone -b sg/readfish-bindings --recursive https://github.com/ratschlab/metagraph.git
cd metagraph
git submodule update --init --recursive
mkdir -p metagraph/build
cd metagraph/build
cmake ..
make -j 16

# setup spumoni
cd $EXPDIR/code
[ -d spumoni ] || git clone --recursive https://github.com/ratschlab/spumoni.git
cd spumoni
mkdir -p build && cd build
cmake ..
make -j 16
make install

# python bindings
cd ..
pip install .

# setup rawhash
cd $EXPDIR/code
[ -d rawhash2 ] || git clone -b cmake --recursive https://github.com/ratschlab/RawHash.git rawhash2
cd rawhash2
git submodule update --init --recursive
mkdir -p build && cd build
cmake NOHDF5=1 NOPOD5=1 ..
make -j 8

# for python bindings
cd ..
pip install .

# setup ReadBouncer
cd $EXPDIR/code
[ -d ReadBouncer ] || git clone --recursive https://github.com/ratschlab/ReadBouncer.git
cd ReadBouncer
git submodule update --init --recursive
mkdir -p build && cd build
cmake ../src
make -j 8

# setup sigmoni
# requires SPUMONI (built above) and Uncalled4
cd $EXPDIR/code
[ -d sigmoni ] || git clone https://github.com/ratschlab/sigmoni.git
cd sigmoni
pip install uncalled4
pip install .

# setup seq2squiggle
cd $EXPDIR/code
[ -d seq2squiggle ] || git clone https://github.com/ZKI-PH-ImageAnalysis/seq2squiggle.git
cd seq2squiggle
pip install .

# Setup Minknow API Simulator
cd $EXPDIR/code
[ -d MinknoApiSimulator ] || git clone https://github.com/ratschlab/MinknoApiSimulator.git
cd MinknoApiSimulator/certs
./generate.sh
cd ..
pip install .

# NOTE: ont-dorado-server is NOT set up by this script. It is ONT's vendored
# Dorado/Guppy basecall server binary distribution, gated behind an ONT
# community login (https://community.nanoporetech.com/). Download the Linux
# server package matching your basecaller version and extract it to
# $EXPDIR/code/ont-dorado-server/ -- see README.md step 5 for checking
# client/server version compatibility.

cd $EXPDIR