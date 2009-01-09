#!/bin/bash

echo INSTALLATION SCRIPT

rm -rf /tmp/pca
mkdir /tmp/pca
export PATH_LLIB=/tmp/pca
export PATH_APLI=/tmp/pca

cp codi/TGZs/fftw-2.1.3.tar.gz $PATH_LLIB
cp codi/TGZs/gnu_licensed_3D_Dock.tar.gz $PATH_APLI
cp codi/TGZs/proteins.tar.gz $PATH_APLI

cd $PATH_LLIB
tar -xvzf fftw-2.1.3.tar.gz
cd fftw-2.1.3
mkdir installation
env CFLAGS="-O3 -march=pentium4" ./configure --prefix=$PWD/installation
make
make install

cd $PATH_APLI
tar -xvzf gnu_licensed_3D_Dock.tar.gz
cd $PATH_APLI/3D_Dock/progs/
tar -zxvf $PATH_APLI/proteins.tar.gz
echo
echo
echo "*****************************************************************"
echo "* Change the FFTW_DIR variable value in the Makefile with that: *"
echo "* FFTW_DIR = $(PATH_LLIB)/fftw-2.1.3/installation               *"
echo "*                                                               *"
echo "* After doing so, just optimize the code and compile using make *"
echo "*****************************************************************"
echo

