#!/bin/bash

export EXPDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $EXPDIR

# Create directory structure
mkdir code logs out results tmp

# setup minimap2
cd $EXPDIR/code
git clone https://github.com/lh3/minimap2.git
cd minimap2
make -j8

# setup collinearity
cd $EXPDIR/code
git clone --recursive https://github.com/ratschlab/collinearity.git
cd collinearity

# for python bindings
pip install .

# to build from source
mkdir build && cd build  
cmake ..  
make -j 8

# setup metagraph
conda install -c bioconda -c conda-forge metagraph
cd $EXPDIR/code
git clone https://github.com/ratschlab/metagraphRF.git
cd metagraphRF
pip install .

# setup spumoni
cd $EXPDIR/code
git clone --recursive https://github.com/ratschlab/spumoni.git
cd spumoni
mkdir build && cd build
cmake ..
make -j 16
make install

# python bindings
cd ..
pip install .

# setup rawhash
cd $EXPDIR/code
git clone -b cmake_merge --recursive https://github.com/ratschlab/RawHash.git rawhash2
cd rawhash2
git submodule update --init --recursive
mkdir -p build && cd build
cmake NOHDF5=1 NOPOD5=1 ..
make -j 8

# for python bindings
cd ..
pip install .

# Setup Minknow API Simulator
cd $EXPDIR/code
git clone https://github.com/ratschlab/MinknoApiSimulator.git
cd MinknoApiSimulator/certs
./generate.sh
cd ..
pip install .

cd $EXPDIR