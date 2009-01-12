#!/bin/bash

echo
echo SCRIPT DE COMPILACION
echo

if [ $# -eq 0 ]; then
   echo "Falta el parametro con el codigo a compilar"
   echo "Los parametros soportados son:"
   echo "     orig       Compila el codigo original"
   echo "     optX       Compila la optimizacion X. X debe ser un valor entre 1 y ...."
else
   if `test $1 = orig`; then
      export BASE_DIR=$PWD
      export PATH_LLIB=$BASE_DIR
      export PATH_APLI=$BASE_DIR
      cd $PATH_APLI
      tar -xvzf $BASE_DIR/TGZs/gnu_licensed_3D_Dock.tar.gz
      cd $PATH_APLI/3D_Dock/progs/
      tar -zxvf $BASE_DIR/TGZs/proteins.tar.gz
      tar -zxvf $BASE_DIR/TGZs/ftdock_outputstests.tar.gz
      rm -rf test?.original
      make clean all
      cd $BASE_DIR
   else
      if `test -d optimizaciones/$1`; then
         cd optimizaciones/$1
         ./compilacion.sh
         cd - >/dev/null
      else
         echo "Codigo a compilar invalido"
         echo "Los parametros soportados son:"
         echo "     orig       Compila el codigo original"
         echo "     optX       Compila la optimizacion X. X debe ser un valor entre 1 y ...."
      fi
   fi
fi
echo

