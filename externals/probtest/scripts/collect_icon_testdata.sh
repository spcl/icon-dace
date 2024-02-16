#!/bin/bash

PROBTEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. &> /dev/null && pwd )"

if [[ -z "$ICON_DIR" ]]; then
    echo "Path to the ICON root directory must be set as 'ICON_DIR'"
    exit 1
fi

if [[ -z "$SUBMIT" ]]; then
    echo "Submit command must be set as 'SUBMIT'"
    exit 1
fi

# generate config file for probtest
# for a global case with multiple output types and a single timestep per output file
EXPNAME=atm_amip_test
python $PROBTEST_DIR/probtest.py init \
    --codebase-install "$ICON_DIR" \
    --experiment-name ${EXPNAME} \
    --file-id NetCDF "*atm_3d*.nc" \
    --file-id NetCDF "*lnd*.nc" \
    --reference probtest_output \
    --member-num 2,5 || exit 1
python $PROBTEST_DIR/probtest.py run-ensemble --submit-command "$SUBMIT" || exit 1

# copy input data
rsync -azP --include="${EXPNAME}*/" --include="*${EXPNAME}*_atm_3d*.nc" --exclude="*" --delete $ICON_DIR/experiments/ icon_data || exit 1
rsync -azP --include="${EXPNAME}*/" --include="*${EXPNAME}*_lnd*.nc" --exclude="*" --delete $ICON_DIR/experiments/ icon_data || exit 1
cp $ICON_DIR/run/exp.${EXPNAME}.run icon_data/${EXPNAME} || exit 1
list_of_logs=( $ICON_DIR/run/LOG.exp.${EXPNAME}.run.*.o )
cp ${list_of_logs[-1]} icon_data/${EXPNAME}/LOG.exp.${EXPNAME}.run.12345678.o || exit 1
cp /users/icontest/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/initial_condition/ifs2icon_1979010100_R02B04_G.nc icon_data/${EXPNAME}/initial_condition.nc || exit 1

# for a LAM case with a single output file that includes all timesteps
EXPNAME=mch_opr_r04b07
python $PROBTEST_DIR/probtest.py init \
    --codebase-install "$ICON_DIR" \
    --experiment-name ${EXPNAME} \
    --file-id NetCDF "*atm_3d_ml*.nc" \
    --file-id NetCDF "*atm_3d_hl*.nc" \
    --file-id NetCDF "*atm_3d_pl*.nc" \
    --reference probtest_output \
    --member-num 2,5 || exit 1
python $PROBTEST_DIR/probtest.py run-ensemble --submit-command "$SUBMIT" || exit 1

# copy input data
rsync -azP --include="${EXPNAME}*/" --include="*${EXPNAME}*_atm_3d*.nc" --exclude="*" --delete $ICON_DIR/experiments/ icon_data || exit 1
cp $ICON_DIR/run/exp.${EXPNAME}.run icon_data/${EXPNAME} || exit 1
cp $ICON_DIR/run/LOG.exp.${EXPNAME}.run.*.o icon_data/${EXPNAME}/LOG.exp.${EXPNAME}.run.12345678.o || exit 1
cp /users/icontest/pool/data/ICON/mch/input/opr_r04b07_lhn_12/ll_2019061512.nc icon_data/${EXPNAME}/initial_condition.nc || exit 1

# copy over the long-running performance benchmark to test visualization
# mkdir -p icon_data/performance
# cp /project/g110/icon-performance-benchmark/daint_gpu_pgi/mch_ch_lowres/* icon_data/performance || exit 1
