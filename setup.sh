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

# metagraph's CLI/server build (BUILD_PYTHON_BINDINGS=OFF) requires Boost with
# static libs. A conda env's boost-cpp package only ships shared libs and its
# CMake config gets found before the system one, so point cmake at the system
# (apt libboost-all-dev) Boost explicitly, installing it first if it's missing.
# (The pymetagraph python bindings below use shared Boost and don't need this.)
BOOST_DIR=$(find /usr/lib/*/cmake -maxdepth 1 -type d -iname 'Boost-*' -print -quit 2>/dev/null)
BOOST_IOSTREAMS_DIR=$(find /usr/lib/*/cmake -maxdepth 1 -type d -iname 'boost_iostreams-*' -print -quit 2>/dev/null)
if [ -z "$BOOST_DIR" ] || [ -z "$BOOST_IOSTREAMS_DIR" ]; then
  echo "metagraph: system Boost (static libs) not found, trying to install libboost-all-dev..."
  if sudo -n true 2>/dev/null; then
    sudo apt-get update && sudo apt-get install -y libboost-all-dev
  elif [ "$(id -u)" = "0" ]; then
    apt-get update && apt-get install -y libboost-all-dev
  fi
  BOOST_DIR=$(find /usr/lib/*/cmake -maxdepth 1 -type d -iname 'Boost-*' -print -quit 2>/dev/null)
  BOOST_IOSTREAMS_DIR=$(find /usr/lib/*/cmake -maxdepth 1 -type d -iname 'boost_iostreams-*' -print -quit 2>/dev/null)
fi

if [ -n "$BOOST_DIR" ] && [ -n "$BOOST_IOSTREAMS_DIR" ]; then
  cmake -DBoost_DIR="$BOOST_DIR" -Dboost_iostreams_DIR="$BOOST_IOSTREAMS_DIR" ..
  make -j 16
else
  echo "WARNING: skipping the metagraph CLI/server build -- could not find or" >&2
  echo "         install a system Boost with static libs (libboost-all-dev)." >&2
  echo "         Ask your admin to run: sudo apt-get install -y libboost-all-dev" >&2
  echo "         then re-run setup.sh to build the metagraph CLI binary." >&2
  echo "         (pymetagraph python bindings will still be installed below.)" >&2
fi

# pymetagraph python bindings (readfish-compatible aligner plugin)
# Separate cmake invocation via scikit-build-core with -DBUILD_PYTHON_BINDINGS=ON,
# which switches metagraph's CMakeLists.txt to shared Boost libs, so the static-lib
# workaround above isn't needed here.
cd $EXPDIR/code/metagraph/metagraph/api/python/pymetagraph
pip install .

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

# pyreadbouncer python bindings
cd $EXPDIR/code/ReadBouncer
pip install .

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