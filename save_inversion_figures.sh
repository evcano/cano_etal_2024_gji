#!/bin/bash
# noise 1: anat 28, fwani 32, fdfwani 31
# noise 2: anat 18, fwani 25, fdfwani 28

ITMAX=28

INVDIR="./inversions_noise_model_2/fdfwani"
OUTDIR=""
TITLE=""

LIM1="550.E+3"
LIM2="350.E+3"

VMIN_1="3000"
VMAX_1="4600"

VMIN_2="-15"
VMAX_2="15"

# DONT EDIT BELOW THIS LINE
# ==================================================
mkdir -p ${OUTDIR}

### data misfit
./plot_data_misfit_evolution.py --fig ${OUTDIR} ${INVDIR}/output/residuals 0 ${ITMAX}

./plot_data_misfit_histograms.py ${INVDIR}/output/residuals 0 --indir2 ${INVDIR}/output/residuals --it2 ${ITMAX} --fig ${OUTDIR}/data_histogram.png

## model misfit
./plot_model_misfit.py --vmin $VMIN_2 --vmax $VMAX_2 --sta ${INVDIR}/specfem2d/DATA/STATIONS --fig ${OUTDIR} inversions_noise_model_1/seisflows_fwd_obs/model_init inversions_noise_model_1/seisflows_fwd_obs/model_init ${INVDIR}/model_init ${INVDIR}/output 0 ${ITMAX} ${LIM1} ${LIM2}

# inverted models
./plot_2dmodel.py --lim1 $LIM1 --lim2 $LIM2 --vmin $VMIN_1 --vmax $VMAX_1 --sta ${INVDIR}/specfem2d/DATA/STATIONS --fig ${OUTDIR}/inverted_00.png  ${INVDIR}/model_init ${INVDIR}/output/MODEL_INIT vs

for ((i=1;i<=${ITMAX};i++))
  do
    j=`printf "%02d" $i`
    ./plot_2dmodel.py --lim1 $LIM1 --lim2 $LIM2 --vmin $VMIN_1 --vmax $VMAX_1 --sta ${INVDIR}/specfem2d/DATA/STATIONS --fig ${OUTDIR}/inverted_${j}.png  ${INVDIR}/model_init ${INVDIR}/output/MODEL_${j} vs
done
