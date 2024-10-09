#!/bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#
# Needs as first line - ksh for DWD, bash for MPI:
#!/bin/bash  #  MPI
#!/bin/ksh   #  DWD
#_____________________________________________________________________________
# Create atmospheric land-sea mask (LSM) for coupled atmosphere-ocean ICON runs:
#
# The LSM from the ocean is strictly logical (0/1) and prescribes the continental
# distribution of the fractional atmospheric LSM. The ocean LSM is interpolated
# with conservative remapping onto the atmospheric grid and added to the extpar
# file as a fraction called cell_sea_land_mask. The logical choices to adopt the
# internal atmosphere LSM (ICON-NWP) to the ocean LSM are defined in ICON within
# mo_ext_data_init.f90.
#
# Info:
# - cell_sea_land_mask (in the ocean grid):
#   long_name = "sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cell"
# - flow chart of LSM in ICON:
#   * lsm_oce2atm_extpar_coupled.sh: interpolate ocean-grid LSM (cell_sea_land_mask)
#                        onto atmo-grid and append to extparfile with same name.
#   * mo_ext_data_init/read_ext_data_atm: read cell_sea_land_mask from extpar file 
#                        and put in variable ext_data%atmo%lsm_ctr_c
#   * mo_ext_data_init/lsm_ocean_atmo: convert with some rules lsm_ctr_c 
#                        to fr_land and fr_lake
#   * mo_atmo_coupling_frame: give ext_data%atm%lsm_ctr_c to YAC representing the atmo LSM
#
# original version for Ruby : Stephan Lorenz and Rene Redler - 2017
# first version for Seamless: Martin Koehler and Rene Redler - 2021-02
# prototyp2              Stephan Lorenz                      - 2022-07-19
# finalize, add to code  Martin Koehler                      - 2022-11-25
# with global ocean grid Stephan Lorenz                      - 2023-01-05
# unify MPI and DWD      Martin Koehler                      - 2023-04-05
# finalize MPI/comments  Stephan Lorenz                      - 2023-04-26
#_____________________________________________________________________________

set -ex

# dwd, levante
site='dwd'
site='levante'

# NWP atmos grids:
atmos_gridID="0030"           # 0012            0030            0024
atmos_refinement="R02B05"     # R02B04          R02B05          R02B06

# ICON-O ocean grids:
ocean_gridID="0035"           # 0036            not             0035          
ocean_refinement="R02B06"     # R02B04          used            R02B06

#_____________________________________________________________________________

case $site in
  'dwd')
     WORKDIR=/hpc/uwork/mkoehler/run-icon/coupled/proto2/lsm-extpar-test
     GRIDDIR=/hpc/rhome/routfox/routfox/icon/grids/public/edzw  # atmospheric grid directory
     EXTPDIR=/hpc/rhome/routfox/routfox/icon/grids/public/edzw  # external parameter directory
     OCENDIR=/hpc/uwork/mkoehler/run-icon/coupled/proto2
     module load cdo/prerelease/2.1.0
     extpname=20180625_tiles       # 20161124_tiles  20180625_tiles  20200917_tiles
     ;;
  'levante')
     WORKDIR=./Create_fractional_lsm_${atmos_gridID}_${ocean_gridID}
     GRIDDIR=/pool/data/ICON/grids/public/edzw                  #  atmospheric NWP grid
     EXTPDIR=/pool/data/ICON/grids/public/edzw                  #  and extpar directory
     OCENDIR=/pool/data/ICON/grids/public/mpim/${ocean_gridID}  #  oceanic MPI grid directory
     set +x
     source /sw/etc/profile.levante
     module load cdo/2.0.6-gcc-11.2.0
     module load nco/5.0.6-gcc-11.2.0
     module list
     set -ex
     extpname=20161124_tiles       # 20161124_tiles available at /pool of MPI
     ;;
esac

#_____________________________________________________________________________

mkdir -p $WORKDIR
cd $WORKDIR

ln -sf ${GRIDDIR}/icon_grid_${atmos_gridID}_${atmos_refinement}_G.nc               atmos_grid.nc
ln -sf ${OCENDIR}/icon_mask_${ocean_gridID}_${ocean_refinement}_G.nc               ocean_mask.nc
ln -sf ${EXTPDIR}/icon_extpar_${atmos_gridID}_${atmos_refinement}_G_${extpname}.nc extpar.nc

# note: ocean_grid only has coastal points and is not suitable for LSM interpolation
#       ocean_mask is global and must be used here for LSM calculation

# not needed if ocean_mask.nc exists:
ln -sf ${OCENDIR}/icon_grid_${ocean_gridID}_${ocean_refinement}_O.nc               ocean_grid.nc
#  create ocean_mask from ocean_grid file by filling missing data into land  - not needed if ocean_mask is global
#cdo setmisstoc,1.0 -selname,cell_sea_land_mask ocean_grid.nc ocean_mask.nc

# first order conservative remapping from ocean to atmospheric grid with fractional LSM [0,1]
# - cdo calculates in 64bit

cdo -b F64 gtc,0                  ocean_mask.nc lsm_0-1.nc           #  ocean mask 0.0/1.0 (gtc,0: >0; 64bit))
cdo -b F64 remapcon,atmos_grid.nc lsm_0-1.nc    lsm_atmos.nc         #  atmos mask 0.0-1.0

# merge new LSM into original extpar file

cdo -O merge  lsm_atmos.nc  extpar.nc                                            lsm_temp.nc
ncatted -O -a rawdata,global,a,c,'GLOBCOVER2009, FAO DSMW, GLOBE, Lake Database' lsm_temp.nc

mv lsm_temp.nc icon_extpar_oceLSM_a${atmos_gridID}_${atmos_refinement}_o${ocean_gridID}_${ocean_refinement}_${extpname}.nc
mv lsm_atmos.nc fractional_atm${atmos_gridID}_oce${ocean_gridID}.nc

\rm -f lsm_0-1.nc

