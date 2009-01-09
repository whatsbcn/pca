#!/bin/bash

echo
echo SCRIPT DE EJECUCION
echo

export BASE_DIR=$PWD
export PATH_LLIB=$BASE_DIR
export PATH_APLI=$BASE_DIR

cd $PATH_APLI/3D_Dock/progs/

echo "./ftdock -static 2pka.parsed -mobile 5pti.parsed > $BASE_DIR/output/test1.out"
./ftdock -static 2pka.parsed -mobile 5pti.parsed > $BASE_DIR/output/test1.out

echo "./ftdock -static 1hba.parsed -mobile 5pti.parsed > $BASE_DIR/output/test2.out"
./ftdock -static 1hba.parsed -mobile 5pti.parsed > $BASE_DIR/output/test2.out

echo "./ftdock -static 4hhb.parsed -mobile 5pti.parsed > $BASE_DIR/output/test3.out"
./ftdock -static 4hhb.parsed -mobile 5pti.parsed > $BASE_DIR/output/test3.out

cd -
