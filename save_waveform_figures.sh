#!/bin/bash

PYCMD="./compare_waveforms.py"

FIGPATH_1=""
OBSDIR_1="inversions_noise_model_1/anat/data_obs"
INIDIR_1="inversions_noise_model_1/anat/output/solver_01"
FINDIR_1="inversions_noise_model_1/anat/output/solver_28"
FLAGS_1="--name --misf --filt 0.02 0.09 --tlim 0 340 --norm --win"

FIGPATH_2=""
OBSDIR_2="inversions_noise_model_1/fwani/data_obs"
INIDIR_2="inversions_noise_model_1/fwani/output/solver_01"
FINDIR_2="inversions_noise_model_1/fwani/output/solver_32"
FLAGS_2="--name --misf --filt 0.02 0.09 --tlim -170 170 --corr --norm --win"

FIGPATH_3=""
OBSDIR_3="inversions_noise_model_1/fdfwani/data_obs"
INIDIR_3="inversions_noise_model_1/fdfwani/output/solver_01"
FINDIR_3="inversions_noise_model_1/fdfwani/output/solver_31"
FLAGS_3="--name --misf --filt 0.02 0.09 --tlim -170 170 --corr --norm --win"

# -------------------------------------------------
# windowing examples
MREC="8P.CAV09"
SREC="8P.CAY12"

# ANAT
FIGN="${FIGPATH_1}/windowing_${SREC}_${MREC}.png"
$PYCMD "${OBSDIR_1}/${MREC}/${SREC}.BXY.semd" "${INIDIR_1}/${MREC}/syn" $FLAGS_1 --fig $FIGN --wind

# FWANI
FIGN="${FIGPATH_2}/windowing_${SREC}_${MREC}.png"
$PYCMD "${OBSDIR_2}/${MREC}/${SREC}.BXY.semd" "${INIDIR_2}/${MREC}/syn" $FLAGS_2 --fig $FIGN --wind

# -------------------------------------------------
# N-S, E-W, NW-SE, NE-SW
MREC_LIST=("8P.CAR02" "8P.CAR07" "8P.CAV09")
SREC_LIST=("8P.CAU03" "8P.CAV07" "8P.CAY12")

for i in {0..2}
do
  MREC=${MREC_LIST[$i]}
  SREC=${SREC_LIST[$i]}

  # ANAT
  FIGN="${FIGPATH_1}/ANAT_${SREC}_${MREC}.png"
  $PYCMD "${OBSDIR_1}/${MREC}/${SREC}.BXY.semd" "${INIDIR_1}/${MREC}/syn" --path3 "${FINDIR_1}/${MREC}/syn" $FLAGS_1 --fig $FIGN
  
  # FWANI
  FIGN="${FIGPATH_2}/FWANI_${SREC}_${MREC}.png"
  $PYCMD "${OBSDIR_2}/${MREC}/${SREC}.BXY.semd" "${INIDIR_2}/${MREC}/syn" --path3 "${FINDIR_2}/${MREC}/syn" $FLAGS_2 --fig $FIGN

  # FD-FWANI
  FIGN="${FIGPATH_3}/FDFWANI_${SREC}_${MREC}.png"
  $PYCMD "${OBSDIR_3}/${MREC}/${SREC}.BXY.semd" "${INIDIR_3}/${MREC}/syn" --path3 "${FINDIR_3}/${MREC}/syn" $FLAGS_3 --fig $FIGN
done
