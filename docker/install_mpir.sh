#!/bin/bash

cd /src

wget mpir.org/mpir-3.0.0.tar.bz2
tar xvjf mpir-3.0.0.tar.bz2
cd mpir-3.0.0

./configure --enable-gmpcompat
make
make install
