#!/bin/ksh
#
# calc_hd_receive_mask.ksh - Calculate hd_receive_mask used in the coupling of HD with ICON
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Trang Van Pham (DWD) and Stefan Hagemann (Hereon)
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# ********* Calculate hd_receive_mask for ICON
# 
# *** Based on the description in Sect. 3.2e in ./docu/readme_coupling.md
#
# 
# ICON grid file, incl. cell_sea_mask:
##ICON_GRID=/pool/data/ICON/edzw-shadow/0030-0035/icon_extpar_oceLSM_a0030_R02B05_o0035_R02B06_20161124_tiles_jsb_cdnc.nc
ICON_GRID=/work/gg0302/g260122/data/icon/trang/grid/icon_grid_a0030_R02B05_o0035_R02B06_20161124_tiles_jsb_cdnc.nc

# Name of output file
OFILE=hd_receive_0030_R02B05.nc

# HD parameter file 
##HD_DIR=/pool/data/ICON/edzw-shadow/indepedent/hd/input/05deg
HD_DIR=/work/gg0302/g260122/HD/input/05deg
HDPARA=${HD_DIR}/hdpara_vs1_11.nc      # HD Parameter File on 0.5 grid - Vs 1.11

# Run directory for script
DRUN=/scratch/g/g260122/tmp

# ************************************
# *** Create the file hd_receive.nc
# ************************************
cd $DRUN

# 1. Select mask with all HD land points
cdo gtc,0. -selvar,FDIR $HDPARA ./mask_hd.nc

# 2. Select ICON points that have land fraction >=0.05
#    In case you want to select all points that have land fraction, 
#    replace 'gec,0.05' by '-gtc,0'
cdo gec,0.05 -selvar,cell_sea_land_mask $ICON_GRID ./iconcoupmask.nc

# 3. Interpolate the ICON coupling land mask with all ICON land/lake points to the HD grid
cdo -L -remapycon,mask_hd.nc ./iconcoupmask.nc ./icon_to_hd_05.nc

# 4. Rename mask variable to hd_receive_mask and store in HD run directory  
cdo setvar,hd_receive_mask icon_to_hd_05.nc ./${OFILE}

rm ./mask_hd.nc ./icon_to_hd_05.nc ./iconcoupmask.nc
