#!/bin/bash

wget https://computation.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz
tar -xzf sundials-2.7.0.tar.gz
cd sundials-2.7.0

mkdir $HOME/sundials-2.7.0/builddir
cd $HOME/sundials-2.7.0/builddir
cmake -DLAPACK_ENABLE=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/sundials-2.7.0 \
      $HOME/sundials-2.7.0

make
make install


ls $HOME/sundials-2.7.0/include