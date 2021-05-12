#!/bin/bash

wget https://github.com/LLNL/sundials/releases/download/v3.1.1/sundials-3.1.1.tar.gz
tar -xzf sundials-3.1.1.tar.gz -C $HOME
cd $HOME/sundials-3.1.1

mkdir $HOME/sundials-3.1.1/builddir
cd $HOME/sundials-3.1.1/builddir
cmake -DLAPACK_ENABLE=ON \
      -DSUNDIALS_INDEX_SIZE=64 \
      -DCMAKE_INSTALL_PREFIX=$HOME/sundials-3.1.1 \
      $HOME/sundials-3.1.1

make
make install


ls $HOME/sundials-3.1.1/include
