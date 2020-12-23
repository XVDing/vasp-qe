#!/bin/bash
root_path=`pwd`
for direction in `ls -F | grep /$`
do
  cd ${root_path}/$direction
     mpirun -np 24 vasp_std
  done
