#!/bin/bash
root_path=`pwd`
for direction in `ls -F | grep /$`
do
  cd ${root_path}/$direction
  for s in strain*
  do
    cd ${root_path}/${direction}$s
       cp ${root_path}/vasp5 ./   # Change the file name of your submit scripts. as me, it is vasp5.
	   sbatch vasp5               # submit your vasp job's command
  done
done
