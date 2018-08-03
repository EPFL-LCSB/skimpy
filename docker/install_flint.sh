#!/bin/bash

cd /src

wget http://www.flintlib.org/flint-2.5.2.tar.gz
tar xzf flint-2.5.2.tar.gz
cd flint-2.5.2

./configure  --with-gmp=/usr/local
make
make install

