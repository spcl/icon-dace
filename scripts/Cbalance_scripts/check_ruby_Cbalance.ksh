#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#-----------------------------------------------------------------------------
# Check carbon conservation in Ruby switch simulations
#
#                                                  Veronika Gayler, Oct. 2020
#-----------------------------------------------------------------------------
# input files needed (in experiment's outdata directory):
#
# ${EXP_ID}_atm_3d_ml_${date}.nc :      mco2vi_phy                  (time avg)
# $(EXP_ID}_lnd_basic_ml_${date}.nc :   carbon_c_sum_natural_ta_box (time avg) 
#                                       carbon_co2_l2a_herb_ta_box  (time avg)
#                                       carbon_co2_l2a_npp_ta_box   (time avg)
#                                       carbon_co2_l2a_resp_ta_box  (time avg)
# $(EXP_ID}_atm_2d_ml_${date}.nc :      co2_flux_lnd                (time avg)
#                                       co2_flux_wtr                (time avg)
# $(EXP_ID}_hamocc_2d_tendencies_${date}.nc  co2flux                (time avg)
# ../bgcout__${date}*Z                 'Global total carbon'        (stepwise) 
#                                      'calcinp'                    (stepwise)
#                                      'orginp'                     (stepwise)
# $(EXP_ID}_cbal_inst_${date}.nc :      carbon_c_sum_natural_ta_box  (inst.)
#                                       mco2vi_phy                   (inst.)
#
# ../bc_land_frac.nc
# ../icon_grid_O.nc                     cell_area
#
#-----------------------------------------------------------------------------
set -x
EXP_ID=vga0300
outdata=/work/mj0060/m220053/ruby-switch/experiments/${EXP_ID}/outdata

cd ${outdata}

YR1=$(ls ${EXP_ID}_lnd_basic_ml_*0101.nc | head -1); YR1=$(echo ${YR1#${EXP_ID}_lnd_basic_ml_} | cut -c-4)
YR2=$(ls ${EXP_ID}_lnd_basic_ml_*0101.nc | tail -1); YR2=$(echo ${YR2#${EXP_ID}_lnd_basic_ml_} | cut -c-4)
DATE="$YR1-$YR2"

if [[ $(cdo showdate ${EXP_ID}_lnd_basic_ml_${YR1}0101.nc | wc -w) == 12 ]]; then
  outfreq=monthly
  muldpm="-muldpm"; tmean="-monmean"; unit="Month"
else
  outfreq=daily
  muldpm=""       ; tmean="-daymean"        ; unit="Day"
fi
#-----------------------------------------------------------------------------
# Preparations
#
# land fraction without lakes and ocean
cdo selvar,land ../bc_land_frac.nc land.nc
# land fraction including lakes
cdo selvar,notsea ../bc_land_frac.nc notsea.nc

# Grid boxes must not have a lake AND an ocean fraction. This might however
# be the case with older grid confgurations. This is a security check, that is
# also done in icon.
cdo selvar,lake ../bc_land_frac.nc lake.nc
cdo ltc,0.9999999999 notsea.nc withsea.nc
cdo mul lake.nc -addc,1 -mulc,-1 withsea.nc lake_corrected.nc
cdo sub notsea.nc lake_corrected.nc land_corrected.nc
if [[ $(cdo -s diff lake.nc lake_corrected.nc) != "" ]]; then
  echo "WARNING: your grid setup is not up-to-date:"
  echo "    Lakes are not allowed in coastal grid cells."
  echo "The land fraction is corrected here in the same way as in the model."
  cp land_corrected.nc land.nc
  cp lake_corrected.nc lake.nc
fi
rm land_corrected.nc lake_corrected.nc

# ocean fraction
cdo selvar,sea ../bc_land_frac.nc sea.nc
# ocean grid cell area [m2]
cdo selvar,cell_area ../icon_grid_O.nc cell_area_oce.nc

#-----------------------------------------------------------------------------
# Atmosphere C inventory
#  mco2vi_phy [kg(CO2) m-2]  (atm_3d_ml)
#
# Convert atmospheric CO2 burden to C content
#          (amco2 = 44.011 in icon ; 44.0095 in echam and jsbach4 comments)
# 44.011 g CO2 = 12.0107 g C  => 1 g CO2 = 12.0107/44.011 g C = 0.272902 g C
# Global surface area: 5.100656e14 m^2
# => 1 kg CO2 / m^2 (area weighted global average) = 5.100656e14 * 0.272902 * 1e-12 Gt C = 139.20 Gt C

cdo -O mergetime -apply,selname,mco2vi_phy "${EXP_ID}_atm_3d_ml_????0101.nc" \
    ${EXP_ID}_atm_3d_ml_mco2vi_phy_${DATE}.nc

cdo -r -setattribute,Catm@units=GtC -setvar,"Catm" -mulc,139.20 -fldmean \
    ${EXP_ID}_atm_3d_ml_mco2vi_phy_${DATE}.nc  ${EXP_ID}_Catm_${DATE}.GtC.nc

#-----------------------------------------------------------------------------
# Land C inventory
#  carbon_c_sum_natural_ta_box [mol(C) m-2]  (lnd_basic_ml)
#
# Convert land C pool from [mol(CO2)/m^2] to [Gt C]
# 1 mol CO2 = 12.0107 g C   => 1 mol CO2 = 1.20107e-14 Gt C
# Global surface area: 5.100656e14 m^2
# => 1 mol CO2 / m^2 (area weighted global average) = 1.20107e-14 * 5.100656e14 Gt C = 6.126245 Gt C

cdo -O mergetime -apply,selname,carbon_c_sum_natural_ta_box "${EXP_ID}_lnd_basic_ml_????0101.nc" \
    ${EXP_ID}_lnd_basic_ml_carbon_c_sum_natural_ta_box_${DATE}.nc

cdo -r -setattribute,Clnd@units=GtC -setvar,"Clnd" -mulc,6.126245 -fldmean -mul notsea.nc  \
    ${EXP_ID}_lnd_basic_ml_carbon_c_sum_natural_ta_box_${DATE}.nc  ${EXP_ID}_Clnd_${DATE}.GtC.nc

#-----------------------------------------------------------------------------
# Ocean C inventory
#  "Global total carbon" [kMol]  (bgcout)
#
#  Convert ocean C content from [kMol C] to [Gt C]
#   1 mol C = 12.0107 g C  =>  1 kMol C = 1.20107e-11 Gt C

# The bgcout file is only written if l_bgc_check=true in hamocc_nml.
# Hamocc prints the inventory checks 4 times per timestep. We are interessted 
# in the last value of the day/month.
#   24 steps per day  x  4 prints per step  => 96 prints per day 
YR=$YR1
rm -f ${EXP_ID}_Coce_${DATE}.tmp ${EXP_ID}_calcinp_${DATE}.tmp ${EXP_ID}_orginp_${DATE}.tmp
while [[ $YR -le $YR2 ]]; do
  grep 'Global total carbon' ../bgcout__${YR}0101T000000Z | cut -f2 -d: >> ${EXP_ID}_Coce_${DATE}.tmp
  grep 'calcinp'             ../bgcout__${YR}0101T000000Z | cut -f2 -d: >> ${EXP_ID}_calcinp_${DATE}.tmp
  grep 'orginp'              ../bgcout__${YR}0101T000000Z | cut -f2 -d: >> ${EXP_ID}_orginp_${DATE}.tmp
 (( YR = YR + 1 ))
done

# Monthly/Daily mean values
sed -n '0~4p' ${EXP_ID}_Coce_${DATE}.tmp > ${EXP_ID}_Coce_${DATE}_hourly.tmp # last print of the timestep (i.e. hour)
cdo -f nc -setattribute,Coce@units=GtC -setvar,"Coce" ${tmean} -settaxis,${YR1}0101,00:00:00,1hour \
    -mulc,1.20107e-11 -input,r1x1 ${EXP_ID}_Coce_${DATE}.GtC.nc < ${EXP_ID}_Coce_${DATE}_hourly.tmp

# Weathering fluxes calc and org [kMol h-1] (i.e. per time step) -> [Gt C] since Jan 01]
#  1 mol(CO2) = 12.0107 gC  =>  1 kMol(CO2) h-1 = 1.20107e-11 GtC h-1
# orginp [kMol(P)], Redfield Ratio Carbon:Phosphate = 122.0

# Calcinp
sed -n '0~4p' ${EXP_ID}_calcinp_${DATE}.tmp > ${EXP_ID}_calcinp_${DATE}_hourly.tmp
cdo -f nc -setattribute,calcinp@units=GtC -setvar,"calcinp" ${tmean} \
    -settaxis,${YR1}0101,00:00:00,1hour -timcumsum -mulc,-1 -mulc,1.20107e-11 \
    -input,r1x1 ${EXP_ID}_calcinp_${DATE}.GtC.nc < ${EXP_ID}_calcinp_${DATE}_hourly.tmp
# Orginp
sed -n '0~4p' ${EXP_ID}_orginp_${DATE}.tmp > ${EXP_ID}_orginp_${DATE}_hourly.tmp
cdo -f nc -setattribute,calcinp@units=GtC -setvar,"orginp" ${tmean} \
    -settaxis,${YR1}0101,00:00:00,1hour -timcumsum -mulc,-1 -mulc,1.20107e-11 -mulc,122 \
    -input,r1x1 ${EXP_ID}_orginp_${DATE}.GtC.nc < ${EXP_ID}_orginp_${DATE}_hourly.tmp


# Instantanous values

# Calcinp
sed -n '0~4p' ${EXP_ID}_calcinp_${DATE}.tmp > ${EXP_ID}_calcinp_${DATE}_hourly.tmp
# get daily sums (and set time steps as in atm output)
cdo -f nc -setattribute,calcinp@units=GtC -setvar,"calcinp" \
    -settime,23:45:00 -shifttime,-1d -seltime,00:00:00 \
    -settaxis,${YR1}0101,01:00:00,1hour -timcumsum -mulc,-1 -mulc,1.20107e-11 \
    -input,r1x1 ${EXP_ID}_calcinp_${DATE}.GtC_daily.nc < ${EXP_ID}_calcinp_${DATE}_hourly.tmp
if [[ ${outfreq} == monthly ]]; then
  # get monthly sums
  cdo -shifttime,-1d -selday,1 -shifttime,1d ${EXP_ID}_calcinp_${DATE}.GtC_daily.nc \
      ${EXP_ID}_calcinp_${DATE}.GtC.nc
else
  mv ${EXP_ID}_calcinp_${DATE}.GtC_daily.nc ${EXP_ID}_calcinp_${DATE}.GtC.nc
fi

# Orginp
sed -n '0~4p' ${EXP_ID}_orginp_${DATE}.tmp > ${EXP_ID}_orginp_${DATE}_hourly.tmp
# get daily sums (and set time steps as in atm output)
cdo -f nc -setattribute,orginp@units=GtC -setvar,"orginp" \
    -settime,23:45:00 -shifttime,-1d -seltime,00:00:00 \
    -settaxis,${YR1}0101,01:00:00,1hour -timcumsum -mulc,-1 -mulc,1.20107e-11 -mulc,122 \
    -input,r1x1 ${EXP_ID}_orginp_${DATE}.GtC_daily.nc < ${EXP_ID}_orginp_${DATE}_hourly.tmp
if [[ ${outfreq} == monthly ]]; then
  # get monthly sums
  cdo -shifttime,-1d -selday,1 -shifttime,1d ${EXP_ID}_orginp_${DATE}.GtC_daily.nc \
      ${EXP_ID}_orginp_${DATE}.GtC.nc
else
  mv ${EXP_ID}_orginp_${DATE}.GtC_daily.nc ${EXP_ID}_orginp_${DATE}.GtC.nc
fi

#-----------------------------------------------------------------------------
# Total carbon
#
# Sum up atmosphere, land and ocean carbon amounts (Ctot is rising due to weathering flux)
#
cdo -O -setattribute,Ctot@units=GtC -setvar,"Ctot" -enssum \
    ${EXP_ID}_Catm_${DATE}.GtC.nc ${EXP_ID}_Clnd_${DATE}.GtC.nc ${EXP_ID}_Coce_${DATE}.GtC.nc \
    ${EXP_ID}_Ctot_${DATE}.GtC.nc

cdo -O -setattribute,Ctot@units=GtC -setvar,"Ctot" -enssum  \
    ${EXP_ID}_Ctot_${DATE}.GtC.nc ${EXP_ID}_calcinp_${DATE}.GtC.nc ${EXP_ID}_orginp_${DATE}.GtC.nc \
    ${EXP_ID}_Ctot-Cinp_${DATE}.GtC.nc

#-----------------------------------------------------------------------------
# C Fluxes
# l2a: co2_flux_lnd [kg m-2 s-1]  (atm_2d_ml)
#      or
#        carbon_co2_l2a_herb_ta_box [kg(CO2) m-2 s-1]  (lnd_basic_ml)
#      + carbon_co2_l2a_npp_ta_box  [kg(CO2) m-2 s-1]  (lnd_basic_ml)
#      + carbon_co2_l2a_resp_ta_box [kg(CO2) m-2 s-1]  (lnd_basic_ml)
#
# Convert [kg(CO2) m-2 s-1] to [GtC d-1]
# 44.0095 g CO2 = 12.0107 g C  => 1 g CO2 = 12.0107/44.0095 g C = 0.272912 g C
# Global surface area: 5.100656e14 m^2
# => 1 kg CO2 / m^2 (area weighted global average) = 5.100656e14 * 0.272912 * 1e-12 Gt C = 139.20 Gt C
#-----------------------------------------------------------------------------

# l2a-flux in ICON-atm output
cdo -O mergetime -apply,selname,co2_flux_lnd "${EXP_ID}_atm_2d_ml_????0101.nc" \
    ${EXP_ID}_atm_2d_ml_co2_flux_lnd_${DATE}.nc
cdo -r -setattribute,Cl2a@units="GtC/${unit}" -setvar,"Cl2a" -mulc,139.20 ${muldpm} -mulc,86400 -fldmean \
    -mul land.nc ${EXP_ID}_atm_2d_ml_co2_flux_lnd_${DATE}.nc  ${EXP_ID}_Cl2a_${DATE}.atm.GtC.nc 

# l2a-flux in ICON-lnd output
cdo -O mergetime -apply,selname,carbon_co2_l2a_herb_ta_box "${EXP_ID}_lnd_basic_ml_????0101.nc" \
    ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_herb_ta_box_${DATE}.nc
cdo -O mergetime -apply,selname,carbon_co2_l2a_npp_ta_box "${EXP_ID}_lnd_basic_ml_????0101.nc" \
    ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_npp_ta_box_${DATE}.nc
cdo -O mergetime -apply,selname,carbon_co2_l2a_resp_ta_box "${EXP_ID}_lnd_basic_ml_????0101.nc" \
    ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_resp_ta_box_${DATE}.nc
cdo -r -setattribute,Cl2a@units="GtC/${unit}" -setvar,"Cl2a" \
    -mulc,139.20 ${muldpm} -mulc,86400 -fldmean -mul notsea.nc \
    -enssum ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_herb_ta_box_${DATE}.nc \
            ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_npp_ta_box_${DATE}.nc  \
            ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_resp_ta_box_${DATE}.nc \
    ${EXP_ID}_Cl2a_${DATE}.lnd.GtC.nc 

# compare the two - should be zero
cdo -r sub ${EXP_ID}_Cl2a_${DATE}.atm.GtC.nc \
           ${EXP_ID}_Cl2a_${DATE}.lnd.GtC.nc \
    test_${EXP_ID}_Cl2a_${DATE}.atm-lnd.GtC.nc


#-----------------------------------------------------------------------------
# 2d-Tests for the land
#-----------------------------------------------------------------------------

# Test 1: C fluxes to the atmosphere: seen by atmosphere vs. seen by land 
#    [kg(CO2) m-2 s-1] -> [kg(CO2) m-2 d-1]

# l2a-flux in ICON-atm output
cdo -setattribute,Cl2a@units="kg(CO2)/m2/${unit}" -setvar,"Cl2a" ${muldpm} -mulc,86400 \
    -mul land.nc ${EXP_ID}_atm_2d_ml_co2_flux_lnd_${DATE}.nc  ${EXP_ID}_Cl2a_${DATE}.atm.nc
cdo -r remapcon,t127grid ${EXP_ID}_Cl2a_${DATE}.atm.nc        ${EXP_ID}_Cl2a_${DATE}.atm.t127.nc

# l2a-flux in ICON-lnd output
cdo -setattribute,Cl2a@units="kg(CO2)/m2/${unit}" -setvar,"Cl2a" ${muldpm} -mulc,86400 \
    -mul notsea.nc -enssum ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_herb_ta_box_${DATE}.nc \
                           ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_npp_ta_box_${DATE}.nc  \
                           ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_resp_ta_box_${DATE}.nc \
    ${EXP_ID}_Cl2a_${DATE}.lnd.nc 

# compare the two - should be zero
cdo remapcon,t127grid ${EXP_ID}_Cl2a_${DATE}.lnd.nc  ${EXP_ID}_Cl2a_${DATE}.lnd.t127.nc

cdo sub ${EXP_ID}_Cl2a_${DATE}.atm.nc ${EXP_ID}_Cl2a_${DATE}.lnd.nc test_${EXP_ID}_Cl2a_atm-lnd.nc
cdo -r remapcon,t127grid test_${EXP_ID}_Cl2a_atm-lnd.nc             test_${EXP_ID}_Cl2a_atm-lnd.t127.nc
rm -f test_${EXP_ID}_Cl2a_atm-lnd.nc


# Test 2: Land C balance: Check cumulative C fluxes to the atmosphere and total land budget 
#
# Total land carbon - instantanous [mol(C) m-2] -> [kg(CO2) m-2]
# 1 mol CO2 = 44.0095 g(CO2) = 0.0440095 kg(CO2)
cdo -O mergetime -apply,selname,carbon_c_sum_natural_ta_box "${EXP_ID}_cbal_inst_ml_????0101.nc" \
    ${EXP_ID}_cbal_inst_ml_carbon_c_sum_natural_ta_box_${DATE}.nc
cdo -setattribute,Clnd_inst@units="kg(CO2)/m2" -setvar,"Clnd_inst" -mulc,0.044011 -mul notsea.nc \
    ${EXP_ID}_cbal_inst_ml_carbon_c_sum_natural_ta_box_${DATE}.nc ${EXP_ID}_Clnd_inst_${DATE}.nc
cdo -r remapcon,t127grid ${EXP_ID}_Clnd_inst_${DATE}.nc           ${EXP_ID}_Clnd_inst_${DATE}.t127.nc

# Cl2a flux [kg(CO2) m-2 d-1] -> Cl2a_sum [kg(CO2) m-2 since Jan 01]
cdo timcumsum ${EXP_ID}_Cl2a_${DATE}.lnd.nc ${EXP_ID}_Cl2a_sum_${DATE}.lnd.nc
cdo -r remapcon,t127grid ${EXP_ID}_Cl2a_sum_${DATE}.lnd.nc ${EXP_ID}_Cl2a_sum_${DATE}.lnd.t127.nc

cdo add ${EXP_ID}_Clnd_inst_${DATE}.nc ${EXP_ID}_Cl2a_sum_${DATE}.lnd.nc \
    test_${EXP_ID}_Clnd_${DATE}.lnd.nc
cdo -r remapcon,t127grid test_${EXP_ID}_Clnd_${DATE}.lnd.nc test_${EXP_ID}_Clnd_${DATE}.lnd.t127.nc
rm -f test_${EXP_ID}_Clnd_${DATE}.lnd.nc

#-----------------------------------------------------------------------------
# Tests for the ocean
#-----------------------------------------------------------------------------

# Test 1: 2D C fluxes to the atmosphere: seen by atmosphere vs. seen by ocean 
#   [kg(CO2) m-2 s-1] -> [kg(CO2) m-2 d-1]

# o2a-flux in ICON-atm output
cdo -O mergetime -apply,selname,co2_flux_wtr "${EXP_ID}_atm_2d_ml_????0101.nc" \
    ${EXP_ID}_atm_2d_ml_co2_flux_wtr_${DATE}.nc
cdo -mulc,86400 -mul sea.nc ${EXP_ID}_atm_2d_ml_co2_flux_wtr_${DATE}.nc \
    ${EXP_ID}_Co2a_${DATE}.atm.nc
cdo -r remapcon,t127grid ${EXP_ID}_Co2a_${DATE}.atm.nc ${EXP_ID}_Co2a_${DATE}.atm.t127.nc

# o2a-flux in ICON-oce output
#   [kMol(C) m-2 s-1] -> [kg(CO2) m-2 d-1]
#   1 kMol(C) = 1 kMol(CO2) = 44.011 kg(CO2)
cdo -O mergetime -apply,selname,co2flux "${EXP_ID}_hamocc_2d_tendencies_????0101.nc" \
    ${EXP_ID}_hamocc_2d_tendencies_co2flux_${DATE}.nc
cdo -mulc,86400 -mulc,44.011 ${EXP_ID}_hamocc_2d_tendencies_co2flux_${DATE}.nc \
    ${EXP_ID}_Co2a_${DATE}.oce.nc 
cdo remapcon,t127grid ${EXP_ID}_Co2a_${DATE}.oce.nc ${EXP_ID}_Co2a_${DATE}.oce.t127.nc

cdo sub ${EXP_ID}_Co2a_${DATE}.atm.t127.nc ${EXP_ID}_Co2a_${DATE}.oce.t127.nc \
    test_${EXP_ID}_Co2a_atm-oce.t127.nc


# Test 2: Ocean C balance: ocean C budget and cumulative C fluxes to and from the ocean
#
# Total ocean carbon - instantanous [kMol(C)] -> [GtC]
#   1 mol C = 12.0107 g C  =>  1 kMol C = 1.20107e-11 Gt C
sed -n '0~96p' ${EXP_ID}_Coce_${DATE}.tmp > ${EXP_ID}_Coce_${DATE}_daily.tmp  # last print of the day (line 96, 192, ...)
cdo -f nc -setattribute,Coce_inst@units=GtC -setvar,"Coce_inst" -settaxis,${YR1}0101,23:45:00,1day \
    -mulc,1.20107e-11 -input,r1x1 ${EXP_ID}_Coce_inst_${DATE}.GtC_daily.nc < ${EXP_ID}_Coce_${DATE}_daily.tmp
if [[ ${outfreq} == monthly ]]; then
  # get the last day of the month
  cdo shifttime,-1d -selday,1 -shifttime,1d ${EXP_ID}_Coce_inst_${DATE}.GtC_daily.nc ${EXP_ID}_Coce_inst_${DATE}.GtC.nc
else
  mv ${EXP_ID}_Coce_inst_${DATE}.GtC_daily.nc ${EXP_ID}_Coce_inst_${DATE}.GtC.nc
fi

# Co2a flux [kMol(C) m-2 s-1] -> Co2a_sum [GtC since Jan 01]
#   1 mol C = 12.0107 g C  =>  1 kMol C = 1.20107e-11 Gt C
cdo timcumsum ${muldpm} -mulc,86400 -mulc,1.20107e-11 -fldsum -mul cell_area_oce.nc \
    ${EXP_ID}_hamocc_2d_tendencies_co2flux_${DATE}.nc  ${EXP_ID}_Co2a_sum_${DATE}.oce.GtC.nc

cdo -O enssum ${EXP_ID}_Coce_inst_${DATE}.GtC.nc ${EXP_ID}_Co2a_sum_${DATE}.oce.GtC.nc \
              ${EXP_ID}_calcinp_${DATE}.GtC.nc   ${EXP_ID}_orginp_${DATE}.GtC.nc \
    test_${EXP_ID}_Coce_${DATE}.GtC.nc


#-----------------------------------------------------------------------------
# clean up
#-----------------------------------------------------------------------------
rm -f *.tmp
rm -f land.nc lake.nc sea.nc withsea.nc notsea.nc cell_area_oce.nc
rm -f ${EXP_ID}_atm_3d_ml_mco2vi_phy_${DATE}.nc
rm -f ${EXP_ID}_atm_2d_ml_co2_flux_wtr_${DATE}.nc
rm -f ${EXP_ID}_atm_2d_ml_co2_flux_lnd_${DATE}.nc
rm -f ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_herb_ta_box_${DATE}.nc
rm -f ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_npp_ta_box_${DATE}.nc
rm -f ${EXP_ID}_lnd_basic_ml_carbon_co2_l2a_resp_ta_box_${DATE}.nc
rm -f ${EXP_ID}_lnd_basic_ml_carbon_c_sum_natural_ta_box_${DATE}.nc
rm -f ${EXP_ID}_cbal_inst_ml_carbon_c_sum_natural_ta_box_${DATE}.nc
rm -f ${EXP_ID}_hamocc_2d_tendencies_co2flux_${DATE}.nc
rm -f ${EXP_ID}_calcinp_${DATE}.GtC_daily.nc ${EXP_ID}_calcinp_${DATE}.GtC.nc
rm -f ${EXP_ID}_orginp_${DATE}.GtC_daily.nc  ${EXP_ID}_orginp_${DATE}.GtC.nc

