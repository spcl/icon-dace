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
