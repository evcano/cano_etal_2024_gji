#!/bin/bash

SPECFEM_DIR=""

date
echo "running directory: `pwd`"
echo


# remove trash
rm -rf bin
rm -rf OUTPUT_FILES

# make directories
mkdir bin
mkdir OUTPUT_FILES

# link executables
cd bin/
ln -s $SPECFEM_DIR/bin/* .
cd ../

# get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# run mesher
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running mesher"
  echo
  ./bin/xmeshfem2D  >> ./OUTPUT_FILES/output_mesher.txt
else
  echo
  echo " running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D >> ./OUTPUT_FILES/output_mesher.txt
fi
if [[ $? -ne 0 ]]; then exit 1; fi

# run simulation 
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running simulation ..."
  echo
  ./bin/xspecfem2D  >> ./OUTPUT_FILES/output_solver.txt
else
  echo
  echo " running simulation on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D >> ./OUTPUT_FILES/output_solver.txt
fi
if [[ $? -ne 0 ]]; then exit 1; fi

rm slurm*.out
echo
date
echo "done"
