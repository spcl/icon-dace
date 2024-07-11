#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#----------------------------------------------------------------------
# JN 11.11.16, script to scale crop and pasture fractions up ("CMIP6 way")
# - scales crop and pasture up with a given veg_ratio_max
# - reduces natural vegetation accordingly -- first reducing primary, then secondary vegetation
#   boundary condition: reduced natural vegetation >= 1%
# -- jsbach does not distinguish primary and secondary vegetation, therefore it is assumed here 
#    that it does not matter which of them is reduced first
# JN 09.03.17 - update: filter with NaNs of the vrm 
# JN 04.09.17 - update: use zero in the vrm as slm 
# (julia.nabel@mpimet.mpg.de)

echo "--> scale_states_with_veg-ratio-max_v2.1-check-0.bash"

inPath=${1}/
outPath=${2}/
inFilePrefix=${3}
outFilePrefix=${4}
vegRatioMaxFile=${5}
vrmVarName=${6}
startYear=${7}
endYear=${8}
calledFrom=${9}
luhRelease=${10}

year=$startYear

while [ $year -le $endYear ]; do
    echo ${year}

    # if year does not have 4 digits...
    if [[ ${year} -le 99  ||  ${year} -gt 9999 ]]; then
      echo "ERROR! Year was ${year} -- this function currently cannot deal with years <100 or >9999"
      exit
    elif [[ ${year} -le 999 ]]; then
      yearString=0${year}
    else
      yearString=${year}
    fi

    #-- test that they initially all sum up to 1 
    cdo -s -add -selvar,gsecd ${inPath}${inFilePrefix}${yearString}.nc -add -selvar,gothr ${inPath}${inFilePrefix}${yearString}.nc -add -selvar,gcrop ${inPath}${inFilePrefix}${yearString}.nc -selvar,gpast ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_init_1.nc
    # TODO: interesting, there are boolean vs float numerical errors -- think through if it might cause problems later

    #-- crate masks for min ratios
    cdo -s -gtc,1e-10 -selvar,${vrmVarName} ${vegRatioMaxFile} ${outPath}TMP_vrm_min-mask.nc
    cdo -s -gtc,1e-10 -selvar,gcrop ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gc_min-mask.nc
    cdo -s -gtc,1e-10 -selvar,gpast ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gp_min-mask.nc
    # - and for 1% nat fraction -> "TR: there should at least be a natural vegetation fraction of 1%
    #   --> where nat fraction is already below 1% we do not scale (i.e. it stays below 1%!)
    cdo -s -gtc,0.01 -add -selvar,gothr ${inPath}${inFilePrefix}${yearString}.nc -selvar,gsecd ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gn-gt-001_mask.nc
    # -> multiply to get common masks where changes are to be done at all
    cdo -s -mul ${outPath}TMP_gn-gt-001_mask.nc ${outPath}TMP_vrm_min-mask.nc ${outPath}TMP_change_mask.nc
    cdo -s -mul ${outPath}TMP_change_mask.nc ${outPath}TMP_gc_min-mask.nc ${outPath}TMP_change_gc_mask.nc
    cdo -s -mul ${outPath}TMP_change_mask.nc ${outPath}TMP_gp_min-mask.nc ${outPath}TMP_change_gp_mask.nc    

    #-- scale up
    cdo -s -div -selvar,gcrop ${inPath}${inFilePrefix}${yearString}.nc -selvar,${vrmVarName} ${vegRatioMaxFile} ${outPath}TMP_gc_scaled.nc
    cdo -s -div -selvar,gpast ${inPath}${inFilePrefix}${yearString}.nc -selvar,${vrmVarName} ${vegRatioMaxFile} ${outPath}TMP_gp_scaled.nc

    #-- only select those that can be changed (else the original unscaled values are used)
    cdo -s -ifthenelse ${outPath}TMP_change_gc_mask.nc ${outPath}TMP_gc_scaled.nc -selvar,gcrop ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gc_changed.nc
    cdo -s -ifthenelse ${outPath}TMP_change_gp_mask.nc ${outPath}TMP_gp_scaled.nc -selvar,gpast ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gp_changed.nc

    #-- TR: there should at least be a natural vegetation fraction of 1%
    cdo -s -lec,0.99 -add ${outPath}TMP_gc_changed.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_sum-le-099_mask.nc
    cdo -s -setname,'crop+pasture_gt_0.99mask' -gtc,0.99 -add ${outPath}TMP_gc_changed.nc ${outPath}TMP_gp_changed.nc ${outPath}${outFilePrefix}${yearString}_scaledCropPlusPasture-gt-099_mask.nc
    cdo -s -subc,0.99 -add ${outPath}TMP_gc_changed.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_sum_diff-to-099.nc
    cdo -s -setname,'crop+pasture_0.99_diff' ${outPath}TMP_sum_diff-to-099.nc ${outPath}${outFilePrefix}${yearString}_scaledCropPlusPasture-diff-to-099.nc

    rm ${outPath}${outFilePrefix}${yearString}_scaledCropPlusPasture*.nc

    # get crop and pasture ratios
    cdo -s -div ${outPath}TMP_gc_changed.nc -add ${outPath}TMP_gc_changed.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_gc_gc-gp_ratio.nc
    cdo -s -div ${outPath}TMP_gp_changed.nc -add ${outPath}TMP_gc_changed.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_gp_gc-gp_ratio.nc

    # if necessary reduce crop and pasture according to their ratio -- but only where changes were done
    cdo -s -sub ${outPath}TMP_gc_changed.nc -mul ${outPath}TMP_gc_gc-gp_ratio.nc ${outPath}TMP_sum_diff-to-099.nc ${outPath}TMP_gc_sub-diff.nc
    cdo -s -ifthenelse ${outPath}TMP_sum-le-099_mask.nc ${outPath}TMP_gc_changed.nc ${outPath}TMP_gc_sub-diff.nc ${outPath}TMP_pre_gc.nc
    cdo -s -ifthenelse ${outPath}TMP_change_gc_mask.nc ${outPath}TMP_pre_gc.nc ${outPath}TMP_gc_changed.nc ${outPath}TMP_gc.nc
    cdo -s -sub ${outPath}TMP_gp_changed.nc -mul ${outPath}TMP_gp_gc-gp_ratio.nc ${outPath}TMP_sum_diff-to-099.nc ${outPath}TMP_gp_sub-diff.nc
    cdo -s -ifthenelse ${outPath}TMP_sum-le-099_mask.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_gp_sub-diff.nc ${outPath}TMP_pre_gp.nc
    cdo -s -ifthenelse ${outPath}TMP_change_gp_mask.nc ${outPath}TMP_pre_gp.nc ${outPath}TMP_gp_changed.nc ${outPath}TMP_gp.nc

    #-- reduce nat veg to fulfil the scaling
    # get the difference between the old and the new sum of crop and pasture, i.e. how much nat veg needs to be reduced
    cdo -s -add ${outPath}TMP_gc.nc ${outPath}TMP_gp.nc ${outPath}TMP_new_sum.nc
    cdo -s -add -selvar,gcrop ${inPath}${inFilePrefix}${yearString}.nc -selvar,gpast ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_old_sum.nc
    cdo -s -sub ${outPath}TMP_new_sum.nc ${outPath}TMP_old_sum.nc ${outPath}TMP_diff_sums.nc

    # reduce primary vegetation (jsbach does not distinguish primary and secondary vegetation, therefore it does not matter which of them is reduced)
    cdo -s -sub -selvar,gothr ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_diff_sums.nc ${outPath}TMP_gothr.nc
    cdo -s -ifthenelse -gec,0 ${outPath}TMP_gothr.nc ${outPath}TMP_gothr.nc -gec,0 ${outPath}TMP_gothr.nc ${outPath}TMP_gothr_gec0.nc
    
    # if primary vegetation was not enought to fulfill scaling reduce secondary
    cdo -s -ifthenelse -ltc,0 ${outPath}TMP_gothr.nc ${outPath}TMP_gothr.nc -ltc,0 ${outPath}TMP_gothr.nc ${outPath}TMP_gothr_lec0.nc
    cdo -s -add -selvar,gsecd ${inPath}${inFilePrefix}${yearString}.nc ${outPath}TMP_gothr_lec0.nc ${outPath}TMP_gsecd.nc

    # TODO: check that they now all still sum up to 1 
    cdo -s -add ${outPath}TMP_gsecd.nc -add ${outPath}TMP_gothr_gec0.nc -add ${outPath}TMP_gc.nc ${outPath}TMP_gp.nc ${outPath}TMP_1_${yearString}.nc

    # merge as new file
    cdo -s merge ${outPath}TMP_gc.nc ${outPath}TMP_gothr_gec0.nc ${outPath}TMP_gp.nc ${outPath}TMP_gsecd.nc ${outPath}TMP_${outFilePrefix}${yearString}.nc

    # - filter with NANs from vrm map
    cdo -s -gec,0 -selvar,${vrmVarName} ${vegRatioMaxFile} ${outPath}TMP_vrm_val-mask.nc
    cdo -s -ifthen ${outPath}TMP_vrm_val-mask.nc ${outPath}TMP_${outFilePrefix}${yearString}.nc ${outPath}TMP_2_${outFilePrefix}${yearString}.nc

    # and with zeros from vrm map
    cdo -s -gtc,0 -selvar,${vrmVarName} ${vegRatioMaxFile} ${outPath}TMP_vrm_non-zero-mask.nc
    cdo -s -ifthen ${outPath}TMP_vrm_non-zero-mask.nc ${outPath}TMP_2_${outFilePrefix}${yearString}.nc ${outPath}${outFilePrefix}${yearString}.nc

    # replace history 
    ncatted -h -O -a '.',global,d,,  ${outPath}${outFilePrefix}${yearString}.nc
    ncatted -h -O -a history,global,o,c,"${thisDate} \n  created with: ${calledFrom} \n from: ${luhRelease}" ${outPath}${outFilePrefix}${yearString}.nc

    (( year = year + 1 ))
done

rm ${outPath}TMP_*
