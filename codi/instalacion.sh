#!/bin/bash

echo INSTALLATION SCRIPT

export BASE_DIR=$PWD
export PATH_LLIB=$BASE_DIR
export PATH_APLI=$BASE_DIR

cd $PATH_LLIB
tar -xvzf $BASE_DIR/TGZs/fftw-2.1.3.tar.gz
cd fftw-2.1.3
mkdir installation
env CFLAGS="-O3 -march=pentium4" ./configure --prefix=$PWD/installation
make
make install

cd $PATH_APLI
tar -xvzf $BASE_DIR/TGZs/gnu_licensed_3D_Dock.tar.gz
cd $PATH_APLI/3D_Dock/progs/
tar -zxvf $BASE_DIR/TGZs/proteins.tar.gz
tar -zxvf $BASE_DIR/TGZs/ftdock_outputstests.tar.gz
make
mkdir $BASE_DIR/output
mv test1.original $BASE_DIR/output/test1.original
mv test2.original $BASE_DIR/output/test2.original
mv test3.original $BASE_DIR/output/test3.original

cd $BASE_DIR

