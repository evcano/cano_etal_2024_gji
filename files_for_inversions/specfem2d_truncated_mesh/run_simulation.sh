#!/bin/bash

SPECFEM_DIR=""

date
echo "running directory: `pwd`"
echo


# remove trash
rm -rf bin
rm -rf OUTPUT_FILES
rm -rf OUTPUT_STEP_1
rm -rf OUTPUT_STEP_2

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

# run noise simulation step 1
./change_noise_step.sh 1

if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running noise simulation step 1..."
  echo
  ./bin/xspecfem2D  >> ./OUTPUT_FILES/output_solver.txt
else
  echo
  echo " running noise simulation step 1 on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D >> ./OUTPUT_FILES/output_solver.txt
fi
if [[ $? -ne 0 ]]; then exit 1; fi

mkdir OUTPUT_STEP_1
mv OUTPUT_FILES/*jpg OUTPUT_STEP_1
mv OUTPUT_FILES/*sem* OUTPUT_STEP_1

# run noise simulation step 2
./change_noise_step.sh 2

if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running noise simulation step 2..."
  echo
  ./bin/xspecfem2D
else
  echo
  echo " running noise simulation step 2 on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D
fi
if [[ $? -ne 0 ]]; then exit 1; fi

mkdir OUTPUT_STEP_2
mv OUTPUT_FILES/*jpg OUTPUT_STEP_2
mv OUTPUT_FILES/*sem* OUTPUT_STEP_2

rm slurm*.out
echo
date
echo "done"
