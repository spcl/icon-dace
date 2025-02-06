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

# -------------------------------------------------------
# Meteogram plots
#
# run as: mtgrm_plot.s 2011010100 R2B04 exp24 dir iFile
#
# attention:
#   scripts need to be local (mtgrm_plot.s, mtgrm_plot.ncl)
#
# Martin Koehler, Dec 2011
# -------------------------------------------------------

#set -ex

dates=${1}
res=${2}
expnum=${3}
echo "mtgrm_plot.s --- Arguments: dates="${dates}" res="${res}" expnum="${expnum}

dir=${4}
iFile=${5}
#dir="/e/uwork/mkoehler/icon/experiments/"${expnum}"/"
#iFile=${dir}"NWP_icon"${res}"_DOM01_"${dates}"_0001_meteogram.nc"
echo "dir = " $dir " iFile = " $iFile

mkdir -p ${dir}"/meteo"
#oType="png" #doesn't work on AIX (NCL 5.2.1)
oType="eps"

set -A iStation  2 3 4 5 6 7 8 9 10 11 12 13 14

set -A varNameSfc \
  P_SFC    PL_Cov   LA_Ind   RO_Dept  Z0       qv_s     w_i      w_snow   \
  TCM      TCH      SHFL     LHFL     RUNOFF_S RUNOFF_G VIO3     HMO3     \
  t_snow   t_s      t_g      FRESHSNW RHO_SNOW H_SNOW   T2M      TD2M     \
  U10M     V10M     SOBT     SOBS     THBS     ALB      RAIN_GSP SNOW_GSP \
  RAIN_CON SNOW_CON 

set -A varName3D \
  P        T        PEXNER   QV       QC       QI       QR       QS       \
  REL_HUM  RHO      THETAV   U        V        CLC      TKVM     TKVH     \
  W        Phalf    t_so     w_so     w_so_ice

set -A iStation   2 
set -A varNameSfc P_SFC    PL_Cov   LA_Ind   RO_Dept  Z0       qv_s     w_i      w_snow   \
  TCM      TCH      SHFL     LHFL     RUNOFF_S RUNOFF_G VIO3     HMO3     \
  t_snow   t_s      t_g      FRESHSNW RHO_SNOW H_SNOW   T2M      TD2M     \
  U10M     V10M     SOBT     SOBS     THBS     ALB      RAIN_GSP SNOW_GSP \
  RAIN_CON SNOW_CON

set -A varName3D    P        T        PEXNER   QV       QC       QI       QR       QS       \
  REL_HUM  RHO      THETAV   U        V        CLC      TKVM     TKVH     \
  W        Phalf    t_so     w_so     w_so_ice

# extended set of variables
# Tobias Goecke (DWD) 02/2018
set -A varName3D          \
P \
T \
PEXNER \
RHO \
THETAV \
U \
V \
W \
TKE \
ddt_tke_hsh \
ddt_tke_pconv \
QV \
QC \
QI \
QR \
QS \
REL_HUM \
QV_DIA \
QC_DIA \
QI_DIA \
CLC \
TKVM \
TKVH \
PHALF \
T_SO \
W_SO \
W_SO_ICE \
W_SO_T_1 \
W_SO_T_2 \
W_SO_T_3 \
W_SO_T_4 \
W_SO_T_5 \
W_SO_T_6 \
T_SO_T_1 \
T_SO_T_2 \
T_SO_T_3 \
T_SO_T_4 \
T_SO_T_5 \
T_SO_T_6


set -A varNameSfc \
PL_COV             \
LA_IND             \
RO_DEPT            \
Z0                 \
QV_S               \
W_I                \
W_SNOW             \
RUNOFF_S           \
RUNOFF_G           \
T_SNOW             \
T_S                \
T_G                \
FRESHSNW           \
RHO_SNOW           \
H_SNOW             \
FR_SEAICE          \
P_SFC              \
TCM                \
TCH                \
SHFL               \
LHFL               \
VIO3               \
HMO3               \
T2M                \
TD2M               \
U10M               \
V10M               \
VBMAX10M           \
dyn_gust           \
con_gust           \
cape_ml            \
SOBT               \
THBT               \
SOBS               \
THBS               \
ALB                \
RAIN_GSP           \
SNOW_GSP           \
RAIN_CON           \
SNOW_CON           \
H_ICE              \
CLCT               \
CLCL               \
CLCM               \
CLCH               \
hbas_con           \
htop_con           \
UMFL_S             \
VMFL_S             \
SWDIFU_S           \
SWDIFD_S           \
PAB_S              \
SWDIR_S            \
T_G_T_1            \
T_G_T_2            \
T_G_T_3            \
T_G_T_4            \
T_G_T_5            \
T_G_T_6            \
T_G_T_7            \
T_G_T_8            \
T_G_T_9            \
SHFL_T_1           \
SHFL_T_2           \
SHFL_T_3           \
SHFL_T_4           \
SHFL_T_5           \
SHFL_T_6           \
SHFL_T_7           \
SHFL_T_8           \
SHFL_T_9           \
LHFL_T_1           \
LHFL_T_2           \
LHFL_T_3           \
LHFL_T_4           \
LHFL_T_5           \
LHFL_T_6           \
LHFL_T_7           \
LHFL_T_8           \
LHFL_T_9           \
SOBS_T_1           \
SOBS_T_2           \
SOBS_T_3           \
SOBS_T_4           \
SOBS_T_5           \
SOBS_T_6           \
SOBS_T_7           \
SOBS_T_8           \
SOBS_T_9           \
THBS_T_1           \
THBS_T_2           \
THBS_T_3           \
THBS_T_4           \
THBS_T_5           \
THBS_T_6           \
THBS_T_7           \
THBS_T_8           \
THBS_T_9           \
FRAC_T_1           \
FRAC_T_2           \
FRAC_T_3           \
FRAC_T_4           \
FRAC_T_5           \
FRAC_T_6           \
FRAC_T_7           \
FRAC_T_8           \
FRAC_T_9           \
snowfrac_t_1       \
snowfrac_t_2       \
snowfrac_t_3       \
snowfrac_t_4       \
snowfrac_t_5       \
snowfrac_t_6       \
TQV                \
TQC                \
TQI                \
TQR                \
TQS                \
TQV_DIA            \
TQC_DIA            \
TQI_DIA            

# -------------------------------------------------------

for station in ${iStation[*]}
do
  for var in ${varNameSfc[*]}
  do
    oFile=${dir}"/meteo/NWP_icon"${res}"_DOM01_"${dates}"_0001_meteogram.loc"${station}"."${var}
    ncl -n mtgrm_plot_sfc.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
      varName=\"${var}\" iStation=${station} expnum=\"${expnum}\"
   #convert -trim -geometry 1000x1000 ${oFile}.pdf ${oFile}.png || true
   #convert -density 100 ${oFile}.eps ${oFile}.png || true
  done	

  for var in ${varName3D[*]}
  do
    oFile=${dir}"/meteo/NWP_icon"${res}"_DOM01_"${dates}"_0001_meteogram.loc"${station}"."${var}
    ncl -n mtgrm_plot.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
      varName=\"${var}\" iStation=${station} expnum=\"${expnum}\"
   #convert -trim -geometry 1000x1000 ${oFile}.pdf ${oFile}.png  || true
   #convert -density 100 ${oFile}.eps ${oFile}.png || true
  done	
done	


exit
