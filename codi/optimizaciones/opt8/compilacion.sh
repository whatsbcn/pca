#!/bin/bash

cp electrostatics.c ../../3D_Dock/progs/electrostatics.c
cp coordinates.c ../../3D_Dock/progs/coordinates.c
cd ../../3D_Dock/progs
make
cd -
