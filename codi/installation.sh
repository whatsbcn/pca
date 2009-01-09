#!/bin/bash

echo INSTALLATION SCRIPT

export BASE_DIR=$PWD
export PATH_LLIB=$BASE_DIR
export PATH_APLI=$BASE_DIR

#cp TGZs/fftw-2.1.3.tar.gz $PATH_LLIB
#cp TGZs/gnu_licensed_3D_Dock.tar.gz $PATH_APLI
#cp TGZs/proteins.tar.gz $PATH_APLI
#cp TGZs/ftdock_outputstests.tar.gz $PATH_APLI

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
mkdir $BASE_DIR/output
mv test1.output $BASE_DIR/outputs/test1.original
mv test2.output $BASE_DIR/outputs/test2.original
mv test3.output $BASE_DIR/outputs/test3.original

cd $BASE_DIR

