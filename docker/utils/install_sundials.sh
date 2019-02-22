#!/bin/bash

#wget https://computation.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz
wget https://computation.llnl.gov/projects/sundials/download/sundials-3.1.1.tar.gz
tar -xzf sundials-3.1.1.tar.gz
cd sundials-3.1.1

mkdir $HOME/sundials-3.1.1/builddir
cd $HOME/sundials-3.1.1/builddir
cmake -DLAPACK_ENABLE=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/sundials-3.1.1 \
      $HOME/sundials-3.1.1

make
make install


ls $HOME/sundials-3.1.1/include