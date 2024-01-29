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

#_____________________________________________________________________________
# Script to construct for a given ICON grid file initial fields and boundary conditions
# for the land model JSBACH4+ and the ICON atmosphere by interpolating existing files
# from ECHAM6 on Gaussian grid
#
# Authors: Reiner Schnur, Veronika Gayler, Rene Redler, MPI-M
#
# Reiner Schnur, MPI-M, 2017-04-10 initial version
#
# Based on scripts from Veronika Gayler, Thomas Raddatz and Marco Giorgetta
#
# The script generates cover fractions relative to the parent tile of the
# respective tiles. This is what ICON-Land generally expects, unless namelist
# parameter 'relative_fractions_in_file = .false.' (jsb_model_nml).
#
# Note: It might be useful to check the jsbach4_ini_files_from_gauss.sh version
#       stored with the intitial files of the previous revision to compare
#       parameter settings.
#
set -eu
#_____________________________________________________________________________

# Settings

coupled=${coupled:-false}         # land sea mask for coupled experiment?
refinement=${refinement:-R02B04}  # ICON-Land grid resolution
bisections=${refinement#*B}  # Determines Gaussian source grid resolution (see below)
atmGridID=${atmGridID:-0043}
oceGridID=${oceGridID:-0035} # Required for coupled configurations
output_dir=${bc_file_dir}_from_gauss  # Output directory for bc ad ic files

start_year=${start_year:-1850}
end_year=${end_year:-1850}

#ggrids="T63"                # additional Gaussian grids for jsbach4 (to run with echam6)
ggrids=""                    # disable Gaussian grids

rootdir=${icon_grid_rootdir:-/pool/data/ICON/grids/public/mpim}
extpar_file=${initial_extpar_file}
[[ ${extpar_file} == "" ]] && \
     extpar_file=/pool/data/ICON/grids/public/mpim/0043/extpar.2016/r0001/icon_extpar_grid_0043_R02B04_G_20210511.nc
icon_grid_file=${icon_grid}
lsm_file=${fractional_mask} # file containing fractional land sea mask
# not used with ICON grids if output_dir is defined
output_root=${output_root_dir:-/pool/data/ICON/grids/public/mpim}
revision=${revision:-r00xx}  # Revision directory that will be generated

debug=false   # true/false
debug_remap_scheme=nn

dryrun=false  # Process files in workdir, but don't copy them to the destination folder (true/false)
clean=true    # Delete workdir after finished (true/false)

# minimum/maximum land grid cell fraction
min_fract=${min_fract:-0.001}
max_fract=${max_fract:-0.999}
# minimum lake fraction
min_lake=0.01

# Currently, only 0 or 11 possible
npfts_list="0 11"

cdo="cdo -s -b F64"

#_____________________________________________________________________________

[[ ${coupled} == "true" ]] && grid_label=$atmGridID-$oceGridID || grid_label=$atmGridID
if [[ ${output_dir} == "" ]]; then
  if [[ ${jsb4icon} == true ]]; then
    output_dir=${output_root}/${grid_label}/land/$revision
  else
    output_dir=${output_root}/${revision}/${grid_name}
  fi
fi

cwddir=$(pwd)
workdir=${work_dir}
[[ $workdir == "" ]] && workdir=${output_dir}/workdir

# Set hd_file to empty if there's no hd_file for this setup available, yet.
hd_file=""

pool_prepare=/pool/data/JSBACH/prepare/                   # pool directory with data needed for initial file generation
grid_path=${rootdir}/$atmGridID

declare -A source_res_dict
source_res_dict[0]=T31
source_res_dict[1]=T31
source_res_dict[2]=T31
source_res_dict[3]=T63
source_res_dict[4]=T63
source_res_dict[5]=T63
source_res_dict[6]=T127
source_res_dict[7]=T127
source_res_dict[8]=T127
source_res_dict[9]=T255
source_res_dict[10]=T255
source_res_dict[11]=T255

#input_root=/pool/data/JSBACH/input

# Note: glac must be at beginning of varlist
varlist="glac lake init_moist snow roughness_length albedo elevation
fao forest_fract maxmoist lai_clim veg_fract surf_temp veg_ratio_max 
albedo_veg_vis albedo_veg_nir albedo_soil_vis albedo_soil_nir roughness_length_oro
soil_depth soil_porosity pore_size_index soil_field_cap heat_capacity heat_conductivity 
moisture_pot hyd_cond_sat wilting_point bclapp fract_org_sl root_depth layer_moist
oromea orostd orosig orogam orothe"
varlist="${varlist} cover_fract"

declare -A sso_dict
sso_dict[oromea]=topography_c
sso_dict[orostd]=SSO_STDH
sso_dict[orosig]=SSO_SIGMA
sso_dict[orogam]=SSO_GAMMA
sso_dict[orothe]=SSO_THETA

#_____________________________________________________________________________

set +u
. ${MODULESHOME}/init/bash
module unload cdo nco netcdf nag
case `hostname` in
  levante*)
    module unload netcdf-c
    module load cdo/2.0.5-gcc-11.2.0
    module load nco/5.0.6-gcc-11.2.0
    module load netcdf-c/4.8.1-gcc-11.2.0
    module load nag/7.1-gcc-11.2.0
    module list
    ;;
   *)
    module load cdo/2.2.0
    module load nco/5.0.1
    module load netcdf-c/4.9.0
    module load nag/7.1
    ;;
esac
set -u

#_____________________________________________________________________________
#
# Function definitions
#_____________________________________________________________________________
#

function generate_gauss_data {

  res_atm=${source_res}
  case ${res_atm} in
    T31) res_oce=GR30
        ;;
    T63) res_oce=GR15
        ;;
    T127) res_oce=TP04
        ;;
    T255) res_oce=TP6M
        ;;
    05) res_oce=""
        ;;
  esac

  if [[ ${jsb4icon} == true ]]; then
    echam_fractional=true
  else
    echam_fractional=false
  fi

  export interactive_ctrl=1                    # 1: swich off interactive mode
  export res_atm=${res_atm}                    # horizontal grid resolution
  export res_oce=${res_oce}                    # ocean model grid (for a coupled setup)
  export ntiles=11                             # number of jsbach tiles
  export dynveg=false                          # setup for simulations with dynamic vegetation
  export c3c4crop=true                         # differentiate between C3 and C4 crops
  export read_pasture=LUH2v2h                  # LUH2v2h / LUH / false
  export pasture_rule=true                     # allocate pastures primarily on grass lands
  export lpasture=true                         # distinguish pastures from grasses
  export year_ct=${year}                       # year the cover_types are derived from (0000 for natural vegetation)
  export year_cf=${year}                       # year cover fractions are derived from (0000 for natural vegetation)
  export landcover_series=false                # generate a series of files with cover_types of different years
  export year_cf2=1859                         # only used with landcover_series
  export echam_fractional=${echam_fractional}  # initial file for echam runs with fractional land sea mask
  export masks_file=default                    # file with land sea mask (default: use echam land sea mask)
  export no_glacier=true                       # ignore glacier mask -> will be later replaced from extpar data
  export pool=/pool/data/ECHAM6/input/r0005/${res_atm}     # directories with echam input data
  if [[ ${res_oce} == dCRUNCEP ]]; then
    export pool=${pool_prepare}/${res_atm}/ECHAM6          # directory with echam input data
  fi
  export pool_land=${pool_prepare}/${res_atm}
  export srcdir=.
  if [[ ${res_atm} == 05 ]]; then
      export read_pasture=LUH
      export pool=/pool/data/JSBACH/prepare/${res_atm}/ECHAM6/
      echo "----------"
      echo " WARNING: LUH2v2h landuse data as used for CMIP6 is not yet available in 0.5 deg resolution."
      echo "----------"
      echo "          We use LUH crop and pasture fractions as in CMIP5, here." 
      echo "          Switch to source grid resolution T127, T63 or T31 if you want LUH2v2h data." 
  fi

  rm -f ../jsbach_*.nc
  # sed -i '/info=/s/false/true/' ./jsbach_init_file.ksh
  ./jsbach_init_file.ksh #vg >& /dev/null
  mv jsbach_*.nc ../

  # Add fields from different source
  if [[ $varlist == *"fract_org_sl"* ]]; then
    ${cdo} remapycon,${source_res}grid ${pool_prepare}05/wise_som-fract_05deg_jsbsoillayers.nc \
        ../gauss_${source_res}_fract_org_sl.nc
  fi

}

function generate_land_sea_mask {

  echo "Generating land sea mask and lake/glacier fractions ..."

  temp=$(mktemp -u --tmpdir=.)

  if [[ ${coupled} == "true" ]]
  then

    # Get fractional mask from fractional mask file depending on ocean grid
    ${cdo} chname,cell_sea_land_mask,notsea -setgrid,${icon_grid_file} \
      -setrtoc,${max_fract},1.1,1 -setrtoc,-1,${min_fract},0 ${lsm_file} ${grid_name}_notsea.nc

    # Glacier mask
    ${cdo} remapycon,${icon_grid_file} gauss_${source_res}_glac.nc ${temp}
    ${cdo} -gtc,0.5   ${temp} ${grid_name}_glac-05_remap.nc
    ${cdo} -gtc,${min_fract} ${temp} ${grid_name}_glac-0001_remap.nc
    # The glacier mask corresponds to the 0.5 mask within the continents and to the "min_fract" mask
    # at the coasts to avoid non-glacier land cells North of Greenland or at the Antarctic coast.
    ${cdo} ifthenelse -gtc,0.7 ${lsm_file}   \
          ${grid_name}_glac-05_remap.nc ${grid_name}_glac-0001_remap.nc \
          ${grid_name}_glac_remap.tmp.nc
    ${cdo} setvar,glac -mul -gtc,0 ${grid_name}_notsea.nc ${grid_name}_glac_remap.tmp.nc ${grid_name}_glac_remap.nc

  else
    # Make sure, the icon grid file includes a valid land sea mask
    if [[ $(cdo -s showvar ${lsm_file} | grep cell_sea_land_mask) == "" ]]; then
      if [[ $(cdo -s showvar ${lsm_file} | grep 'FR_LAND') == "" ]]; then
        echo "Land sea mask not found: lsm_file ${lsm_file} does not contain 'cell_sea_land_mask' nor 'FR_LAND'."
        exit 1
      else
        # Get land sea and glacier masks from initial extpar file
        if [[ $(cdo -s showvar ${lsm_file} | grep 'FR_LAND ') != "" ]]; then
          # Extpar file contains FR_LAND, FR_LAKE and ICE (typically in NWP extpar file)
          lsm_var=FR_LAND
          ${cdo} setname,notsea -setrtoc,-1,${min_fract},0 -setrtoc,${max_fract},1.1,1 -setmisstoc,0 -setmissval,-9.e33 \
             -setgrid,${icon_grid_file} -add -selvar,FR_LAND ${extpar_file} -selvar,FR_LAKE ${extpar_file} \
             ${grid_name}_notsea.nc
          ${cdo} setname,glac -setgrid,${icon_grid_file} -mul -gtc,0 ${grid_name}_notsea.nc \
             -selvar,ICE ${extpar_file} ${grid_name}_glac_remap.nc
        elif [[ $(cdo -s showvar ${lsm_file} | grep 'FR_LAND_TOPO') != "" ]]; then
          # Extpar file contains FR_LAND_TOPO including lake fractions (typically in MPIM extpar file)
          lsm_var=FR_LAND_TOPO
          ${cdo} setname,notsea -setrtoc,-1,${min_fract},0 -setrtoc,${max_fract},1.1,1 -setmisstoc,0 -setmissval,-9.e33 \
             -setgrid,${icon_grid_file} -selvar,FR_LAND_TOPO ${extpar_file} ${grid_name}_notsea.nc
          # Just a dummy glacier mask: will be replaced by extpar mask ...
          ${cdo} mulc,0. ${grid_name}_notsea.nc ${grid_name}_glac_remap.nc
        else
          echo "Land sea mask in ${extpar_file} not found. It is neither FR_LAND nor FR_LAND_TOPO"
          exit 1
        fi
      fi
    else
      # Get land sea masks from grids file (This results in non-fractional 0/1 land sea mask!)
      lsm_var=cell_sea_land_mask
      ${cdo} chname,cell_sea_land_mask,notsea -setrtoc,${max_fract},2.1,1 -setrtoc,-2,${min_fract},0 \
          -selvar,cell_sea_land_mask ${lsm_file} ${grid_name}_notsea.nc
      # Just a dummy glacier mask: will be replaced by extpar mask ...
      ${cdo} mulc,0. ${grid_name}_notsea.nc ${grid_name}_glac_remap.nc
    fi
    if [[ $(cdo -s output -fldsum -gtc,0 -selvar,${lsm_var} ${lsm_file} | tr -d ' ') == 0 ]]
    then
      echo "Invalid land sea mask in lsm_file ${lsm_file}:"
      echo "Variable cell_sea_land_mask does not have any land grid cells."
      echo "Check if there is another grid file, e.g. with suffix _Glsm.nc."
      exit 1
    fi
  fi
  ncatted -a long_name,notsea,m,c,'Fraction of land+lake' ${grid_name}_notsea.nc

  # Generate lake fractions, making sure that fractional coastal grid points don't contain any lakes
  # Set small lake fractions below "min_fract" to zero
  ${cdo} setvar,lake -setrtoc,-1,${min_lake},0 gauss_${source_res}_lake.nc gauss_${source_res}_lake_ext.nc
  ${cdo} remapdis,${grid_name}_notsea.nc gauss_${source_res}_lake_ext.nc ${temp}
  ${cdo} setvar,lake -setmissval,-9.e33 -ifthen -ltc,0.5 ${grid_name}_glac_remap.nc -ifthen -gec,${max_fract} ${grid_name}_notsea.nc \
	 -setrtoc,-1,${min_lake},0 ${temp} ${grid_name}_lake_remap.nc

  rm ${temp}

  # Used to mask out ocean with missing values, debug_file needs this
  ${cdo} setctomiss,0 -setmissval,-9e33 ${grid_name}_notsea.nc ${grid_name}_notsea_miss.nc
  debug_file ${grid_name}_notsea_miss.nc
  debug_file ${grid_name}_notsea.nc
  debug_file ${grid_name}_lake_remap.nc
  debug_file ${grid_name}_glac_remap.nc
}

function extrapolate_var {

  local var=$1     # Variable name
  local cdo_files="gauss_${source_res}_${var}.nc gauss_${source_res}_${var}_ext.nc"
  if [[ -f gauss_${source_res}_${var}.nc ]]
  then
    temp=$(mktemp -u --tmpdir=.)
    ${cdo} setmissval,-9.e33 gauss_${source_res}_${var}.nc ${temp}
    mv ${temp} gauss_${source_res}_${var}.nc
  fi
  case ${var} in
    init_moist | layer_moist | roughness_length | root_depth | fao )
      echo "Extrapolating $var ..."
      ${cdo} setmisstoc,0 -setmissval,-9.e33 -fillmiss -setmissval,0 ${cdo_files}
      ;;
    albedo )
      echo "Extrapolating $var ..."
      ${cdo}  setmissval,-9.e33 -fillmiss -setmissval,0.07 ${cdo_files}
      ;;
    albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
      echo "Extrapolating $var ..."
      ${cdo} copy ${cdo_files}
      ;;
    maxmoist )
      echo "Extrapolating $var ..."
      ${cdo}  setmissval,-9.e33 -fillmiss -setmissval,0 ${cdo_files}
      ;;
    forest_fract| snow )
      echo "Extrapolating $var ..."
      ${cdo} setvar,${var} -fillmiss -ifthen gauss_${source_res}_notsea_miss.nc ${cdo_files}
      ;;
    lai_clim | veg_fract | veg_ratio_max | roughness_length_oro | \
    soil_depth | bclapp | soil_field_cap | heat_capacity | heat_conductivity | hyd_cond_sat | \
    moisture_pot | pore_size_index | soil_porosity | wilting_point | fract_org_sl )
      echo "Extrapolating $var ..."
      # Do not use values from glacier grid cells for extrapolation
      #  (fillmiss does not fill all missing values, that's why setmisstoc is also needed.)
      ${cdo} setvar,${var} -setmisstoc,0 -fillmiss -ifthen gauss_${source_res}_non-glac-land.nc \
      ${cdo_files}
      ;;
    cover_fract )
      echo "Extrapolating $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} splitlevel gauss_${source_res}_${var}.nc cfract
      ${cdo} sub cfract000001.nc gauss_${source_res}_glac.nc cfract_tmp.nc && mv cfract_tmp.nc cfract000001.nc
      ${cdo} merge cfract0000??.nc ${temp} && rm cfract0000??.nc
      ${cdo} setvar,${var} -setmisstoc,0 -setmissval,-9.e33 -fillmiss -ifnotthen gauss_${source_res}_glac.nc \
             -ifthen gauss_${source_res}_notsea_miss.nc  ${temp} gauss_${source_res}_${var}_ext.nc
      rm ${temp}
      ;;
    #* ) echo "Not extrapolating $var ..."
    #   ;;
  esac
  return 0

}

function remap_var {

  local var=$1
  case ${var} in
    fao )
      echo "Remapping $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remapnn,${icon_grid_file} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setname,fao -setmisstoc,0 -setmissval,-9e33 -ifthen ${grid_name}_notsea_miss.nc \
             ${temp} ${grid_name}_${var}_remap.nc
      rm ${temp}
      ;;
    albedo )
      echo "Remapping $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setvar,${var} -setmisstoc,0.07 -ifthenelse ${grid_name}_glac_remap.nc \
             -mulc,0.7 ${grid_name}_glac_remap.nc ${temp} ${grid_name}_${var}_remap.nc
      rm ${temp}
      ;;
    albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir )
      echo "Remapping $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setvar,${var} ${temp} ${grid_name}_${var}_remap.nc
      rm ${temp}
      ;;
    elevation )
      echo "Extracting elevation from EXTPAR file ..."
      ${cdo} setvar,elevation \
            -setgrid,${icon_grid_file} -selvar,topography_c ${extpar_file} \
            ${grid_name}_${var}_remap.nc
      ;;
    oro* )
      echo "Extracting $var from EXTPAR file ..."
      # Note: setrtoc in followoing is necessary because SSO variables in extpar file contain missing values
      ${cdo} setvar,${var} \
            -setgrid,${icon_grid_file} -setrtoc,-9e40,-1000,0 -setrtoc,100000,9e40,0 \
            -selvar,${sso_dict[${var}]} ${extpar_file} \
            ${grid_name}_${var}_remap.nc
      ;;
    cover_fract )
      # npft should be consistent with cover_fract from the JSBACH3 input file
      echo "Remapping $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
      ${cdo} splitlevel ${grid_name}_${var}_remap.nc cfract
      local n=1
      #while [[ ${n} -le ${npfts} ]]
      while [[ ${n} -le 11 ]]
      do
        nt=$(printf "%2.2i" $n)
        ${cdo} setmissval,-9.e33 -mul cfract0000${nt}.nc -mulc,-1 -subc,1 ${grid_name}_glac_remap.nc cfract${nt}.nc && rm cfract0000${nt}.nc
        (( n = n + 1 ))
      done
      ${cdo} setmissval,-9.e33 -add cfract01.nc ${grid_name}_glac_remap.nc ${temp} && mv ${temp} cfract01.nc
      rm ${grid_name}_${var}_remap.nc
      ${cdo} merge cfract??.nc ${grid_name}_${var}_remap.nc
      ;;
    roughness_length | roughness_length_oro )
      echo "Remapping $var ..."
      [[ "${var}" == "roughness_length" ]] && oce_rough=0.001 || oce_rough=0.0
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setvar,${var} -setmissval,-9.e33 -setmisstoc,${oce_rough} -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
      rm ${temp}
      ;;
    glac | lake )
     # In contreast to case *) glacier values are not set to zero. This assures somehow meaningful values
     # in case the extpar glacier mask is different from the glacier mask used here.
      ;;
    surf_temp )
      echo "Remapping $var ..."
      temp=$(mktemp -u --tmpdir=.)
      ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
      ${cdo} setvar,${var} ${temp} ${grid_name}_${var}_remap.nc
      rm ${temp}
      ;;
   soil_depth | bclapp | soil_field_cap | heat_capacity | heat_conductivity | hyd_cond_sat | \
   moisture_pot | pore_size_index | soil_porosity | wilting_point | fract_org_sl )
     echo "Remapping $var ..."
     temp=$(mktemp -u --tmpdir=.)
     ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
     case ${var} in
       soil_depth )
         ${cdo} setvar,${var} -setmisstoc,0.5 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       bclapp )
         ${cdo} setvar,${var} -setmisstoc,4.5 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       soil_field_cap )
         ${cdo} setvar,${var} -setmisstoc,0.229 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       heat_capacity )
         ${cdo} setvar,${var} -setmisstoc,2.e+6 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       heat_conductivity )
         ${cdo} setvar,${var} -setmisstoc,7. -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       hyd_cond_sat )
         ${cdo} setvar,${var} -setmisstoc,5.e-6 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       moisture_pot )
         ${cdo} setvar,${var} -setmisstoc,0.15 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       pore_size_index )
         ${cdo} setvar,${var} -setmisstoc,0.2 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       soil_porosity )
         ${cdo} setvar,${var} -setmisstoc,0.45 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       wilting_point )
         ${cdo} setvar,${var} -setmisstoc,0.15 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       fract_org_sl )
         ${cdo} setvar,${var} -setmisstoc,0.5 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
       * )
         ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
         ;;
     esac
     rm ${temp}
     ;;
    * )
      if [[ -f gauss_${source_res}_${var}_ext.nc ]]
      then
        echo "Remapping $var ..."
        temp=$(mktemp -u --tmpdir=.)
        ${cdo} remap,${icon_grid_file},${remap_weights} gauss_${source_res}_${var}_ext.nc ${temp}
        ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea_miss.nc ${temp} ${grid_name}_${var}_remap.nc
        ${cdo} mul ${grid_name}_${var}_remap.nc -mulc,-1 -subc,1 ${grid_name}_glac_remap.nc ${temp}
        mv ${temp} ${grid_name}_${var}_remap.nc
        #rm ${temp}
      fi
      ;;
  esac
  # if [[ ${var} == "glac" ]]
  # then
  #     temp=$(mktemp -u --tmpdir=.)
  #     ${cdo} gec,0.5 ${grid_name}_glac_remap.nc ${temp}
  #     mv ${temp} ${grid_name}_glac_remap.nc
  # fi
  # Copy attributes but delete CDI_grid_type otherwise it messes up icon grid files (only necessary with new nco version on levante)
  if [ -f gauss_${source_res}_${var}.nc ]; then
      ncks -A -C -H -v .${var} gauss_${source_res}_${var}.nc ${grid_name}_${var}_remap.nc >& /dev/null
      ncatted -O -a CDI_grid_type,,d,, ${grid_name}_${var}_remap.nc >& /dev/null
  fi
  [[ -f ${grid_name}_${var}_remap.nc ]] && debug_file ${grid_name}_${var}_remap.nc
  return 0
}

function pseudo_remap_var {

  local var=$1     # Variable name

  [[ ! -f ${grid_name}_notsea.nc ]] && ${cdo} setvar,notsea gauss_${grid_name}_slm.nc ${grid_name}_notsea.nc
  if [[ ${var} == oro* ]]; then
  	VAR=$(echo ${var} | tr -s 'a-z' 'A-Z')
	  ${cdo} -setvar,${var} -selvar,${VAR}  ${pool}/${res_atm}${res_oce}_jan_surf.nc  ${grid_name}_${var}_remap.nc
  else
	  ln -s gauss_${grid_name}_${var}.nc ${grid_name}_${var}_remap.nc
  fi

}

function generate_fractions {

  ${cdo} setvar,fract_glac -setmisstoc,0 -gec,0.5 ${grid_name}_glac_remap.nc ${grid_name}_fract_glac.nc
  ncatted -a long_name,fract_glac,m,c,'Fraction of glacier tile rel. to land tile' ${grid_name}_fract_glac.nc
  debug_file ${grid_name}_fract_glac.nc
  ${cdo} setvar,fract_veg -ltc,0.5 ${grid_name}_fract_glac.nc ${grid_name}_fract_veg.nc
  ncatted -a long_name,fract_veg,m,c,'Fraction of vegetated tile rel. to land tile' ${grid_name}_fract_veg.nc
  debug_file ${grid_name}_fract_veg.nc
  ${cdo} setvar,sea -setrtoc,-1,0,0 -mulc,-1 -subc,1 ${grid_name}_notsea.nc ${grid_name}_sea.nc
  ncatted -a long_name,sea,m,c,'Fraction of ocean' ${grid_name}_sea.nc
  debug_file ${grid_name}_sea.nc
  ${cdo} setvar,fract_lake -setmisstoc,0 ${grid_name}_lake_remap.nc ${grid_name}_fract_lake.nc
  ncatted -a long_name,fract_lake,m,c,'Fraction of lake tile rel. to box tile' ${grid_name}_fract_lake.nc
  debug_file ${grid_name}_fract_lake.nc
  ${cdo} setvar,fract_land -mulc,-1 -subc,1 ${grid_name}_fract_lake.nc ${grid_name}_fract_land.nc
  ncatted -a long_name,fract_land,m,c,'Fraction of land tile rel. to box tile' ${grid_name}_fract_land.nc
  debug_file ${grid_name}_fract_land.nc
  # The following three are duplicates for compatibility with echam_phy_init, should be removed at a later time.
  ${cdo} setvar,glac ${grid_name}_fract_glac.nc ${grid_name}_glac.nc
  ncatted -a long_name,glac,m,c,'Fraction of glacier' ${grid_name}_glac.nc
  debug_file ${grid_name}_glac.nc
  ${cdo} setvar,lake -mul ${grid_name}_notsea.nc ${grid_name}_fract_lake.nc ${grid_name}_lake.nc
  ncatted -a long_name,lake,m,c,'Fraction of lakes' ${grid_name}_lake.nc
  debug_file ${grid_name}_lake.nc
  ${cdo} setvar,land -mul ${grid_name}_notsea.nc ${grid_name}_fract_land.nc ${grid_name}_land.nc
  ncatted -a long_name,land,m,c,'Fraction of land' ${grid_name}_land.nc
  debug_file ${grid_name}_land.nc

  # Generate PFT fractions
  local i=1
  for npfts in $npfts_list
  do
    if [[ $npfts -gt 0 ]]
    then
      pft_tag="${npfts}pfts_"
      i=1
      while [[ ${i} -le ${npfts} ]]
      do
        ipft=$(printf "%2.2i" $i)
        ${cdo} chname,cover_fract,pft${ipft} -sellevel,$i ${grid_name}_cover_fract_remap.nc pft${ipft}.nc
        if [[ ${i} -eq 1 ]]
        then
          ${cdo} setrtoc,1,1.1,0 pft${ipft}.nc zpft && mv zpft pft${ipft}.nc
        fi
        ncwa -a ntiles pft${ipft}.nc zpft
        ${cdo} setmisstoc,0 zpft pft${ipft}.nc && rm zpft
        (( i = i + 1 ))
      done
      ${cdo} enssum pft*.nc fract_sum.nc
      ${cdo} mulc,-1 -subc,1 -gtc,0 -setmisstoc,0 fract_sum.nc fract_sum2.nc && mv fract_sum2.nc fract_sum.nc
      ${cdo} add pft03.nc fract_sum.nc zpft && mv zpft pft03.nc
      rm fract_sum.nc
      for i in $(seq 1 $npfts)
      do
        ipft=$(printf "%2.2i" $i)
        ${cdo} setname,fract_pft${ipft} pft${ipft}.nc ${grid_name}_fract_pft${ipft}_remap.nc && rm pft${ipft}.nc
        ncatted -a long_name,fract_pft${ipft},m,c,"Fraction of pft${ipft} tile rel. to veg tile" ${grid_name}_fract_pft${ipft}_remap.nc
        debug_file ${grid_name}_fract_pft${ipft}_remap.nc
      done
    fi
  done
}

function crosscheck_vars {

  # Set minimum soil depth of non-glacier grid cells to 0.01
  temp=$(mktemp -u --tmpdir=.)
  ${cdo} ifthenelse -mulc,-1 -subc,1 ${grid_name}_fract_glac.nc -setrtoc,0,0.1,0.1 ${grid_name}_soil_depth_remap.nc ${grid_name}_soil_depth_remap.nc ${temp}
  mv ${temp} ${grid_name}_soil_depth_remap.nc
  debug_file ${grid_name}_soil_depth_remap.nc

  # Calculate maximum soil water amount from soil depth and field capacity
  temp=$(mktemp -u --tmpdir=.)
  ${cdo} mul ${grid_name}_soil_field_cap_remap.nc ${grid_name}_soil_depth_remap.nc max_soilwater.nc
  ${cdo} min ${grid_name}_maxmoist_remap.nc max_soilwater.nc ${temp}
  mv ${temp} ${grid_name}_maxmoist_remap.nc
  debug_file ${grid_name}_maxmoist_remap.nc
  rm max_soilwater.nc

  # Make sure initial soil moisture is less then maxmoist
  temp=$(mktemp -u --tmpdir=.)
  ${cdo} -le ${grid_name}_init_moist_remap.nc ${grid_name}_maxmoist_remap.nc le-mask.tmp
  ${cdo} ifthenelse le-mask.tmp ${grid_name}_init_moist_remap.nc ${grid_name}_maxmoist_remap.nc ${temp}
  mv ${temp} ${grid_name}_init_moist_remap.nc

  # Set missing value of surface temperature to 280K
  temp=$(mktemp -u --tmpdir=.)
  ${cdo} setmisstoc,280 ${grid_name}_surf_temp_remap.nc ${temp}
  mv ${temp} ${grid_name}_surf_temp_remap.nc
  debug_file ${grid_name}_surf_temp_remap.nc
}

function copy_uuid {

  if [[ ${jsb4icon} == true ]]; then
    local file=$1

    uuidOfHGrid=$(ncdump -h ${icon_grid_file} | grep ':uuidOfHGrid = ' | cut -f2 -d'"')
    ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} ${file}
  fi

}

function generate_output_files {

  local outfile=""
  local pft_tag=""

  temp=$(mktemp -u --tmpdir=.)

  for npfts in $npfts_list
  do
    [[ $npfts -gt 0 ]] && pft_tag="${npfts}pfts_" || pft_tag=""
    outfile="${grid_name}_bc_land_frac_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    if [[ ${npfts} -gt 0 ]]
    then
      ${cdo} merge \
          ${grid_name}_notsea.nc \
          ${grid_name}_fract_glac.nc \
          ${grid_name}_sea.nc \
          ${grid_name}_fract_lake.nc \
          ${grid_name}_fract_land.nc \
          ${grid_name}_fract_veg.nc \
          ${grid_name}_veg_ratio_max_remap.nc \
          ${grid_name}_fract_pft*_remap.nc \
          ${grid_name}_land.nc \
          ${grid_name}_lake.nc \
          ${grid_name}_glac.nc \
          ${outfile}
    else
      ${cdo} merge \
          ${grid_name}_notsea.nc \
          ${grid_name}_fract_glac.nc \
          ${grid_name}_sea.nc \
          ${grid_name}_fract_lake.nc \
          ${grid_name}_fract_land.nc \
          ${grid_name}_fract_veg.nc \
          ${grid_name}_veg_ratio_max_remap.nc \
          ${grid_name}_land.nc \
          ${grid_name}_lake.nc \
          ${grid_name}_glac.nc \
          ${outfile}
    fi
    comment=$(cdo -s showattribute,comment ${input_file} | sed -e 's/\\n//g' -e 's/\"//g' | tr '\n' '; ' | tr -s " " | cut -d "=" -f2)
    ${cdo} setattribute,comment="${comment}" ${outfile} ${temp} && mv ${temp} ${outfile}

    #ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

  done

  outfile="${grid_name}_bc_land_phys_${year}.nc"
  echo "Generating output file ${outfile} ..."
  ${cdo} merge \
      ${grid_name}_lai_clim_remap.nc \
      ${grid_name}_veg_fract_remap.nc \
      ${grid_name}_roughness_length_remap.nc \
      ${grid_name}_roughness_length_oro_remap.nc \
      ${grid_name}_albedo_remap.nc \
      ${grid_name}_albedo_veg_vis_remap.nc \
      ${grid_name}_albedo_veg_nir_remap.nc \
      ${grid_name}_albedo_soil_vis_remap.nc \
      ${grid_name}_albedo_soil_nir_remap.nc \
      ${grid_name}_forest_fract_remap.nc \
      ${outfile}
  #ncrename -d cell,ncells ${outfile} >& /dev/null
  copy_uuid ${outfile}
  debug_file ${outfile}

  outfile="${grid_name}_bc_land_soil_${year}.nc"
  echo "Generating output file ${outfile} ..."

  # Initial files on the Gaussian grid should not contain missing values. Variable
  # fract_org_sl is adapted here (after extrapolation), as the icon grid initial
  # files should not be affected.
  if [[ $(echo ${grid_name} | cut -c1) == T ]]; then   # Gausian grid
    cdo setmisstoc,0 ${grid_name}_fract_org_sl_remap.nc  ${grid_name}_fract_org_sl_remap.tmp
    mv ${grid_name}_fract_org_sl_remap.tmp ${grid_name}_fract_org_sl_remap.nc
  fi

  ${cdo} merge \
      ${grid_name}_soil_depth_remap.nc \
      ${grid_name}_root_depth_remap.nc \
      ${grid_name}_fract_org_sl_remap.nc \
      ${grid_name}_bclapp_remap.nc \
      ${grid_name}_soil_field_cap_remap.nc \
      ${grid_name}_heat_capacity_remap.nc \
      ${grid_name}_heat_conductivity_remap.nc \
      ${grid_name}_hyd_cond_sat_remap.nc \
      ${grid_name}_moisture_pot_remap.nc \
      ${grid_name}_pore_size_index_remap.nc \
      ${grid_name}_soil_porosity_remap.nc \
      ${grid_name}_wilting_point_remap.nc \
      ${grid_name}_fao_remap.nc \
      ${grid_name}_maxmoist_remap.nc \
      ${outfile}
  #ncrename -d cell,ncells ${outfile} >& /dev/null
  copy_uuid ${outfile}
  debug_file ${outfile}

  outfile="${grid_name}_bc_land_sso_${year}.nc"
  echo "Generating output file ${outfile} ..."
  ${cdo} merge \
      ${grid_name}_elevation_remap.nc \
      ${grid_name}_oromea_remap.nc \
      ${grid_name}_orostd_remap.nc \
      ${grid_name}_orosig_remap.nc \
      ${grid_name}_orogam_remap.nc \
      ${grid_name}_orothe_remap.nc \
      ${outfile}
  #ncrename -d cell,ncells ${outfile} >& /dev/null
  copy_uuid ${outfile}
  debug_file ${outfile}

  outfile="${grid_name}_ic_land_soil_${year}.nc"
  echo "Generating output file ${outfile} ..."
  ${cdo} merge \
      ${grid_name}_init_moist_remap.nc \
      ${grid_name}_surf_temp_remap.nc \
      ${grid_name}_layer_moist_remap.nc \
      ${grid_name}_snow_remap.nc \
      ${outfile}
  #ncrename -d cell,ncells ${outfile} >& /dev/null
  copy_uuid ${outfile}
  debug_file ${outfile}

}

function generate_frac_output_files {

  local outfile=""
  local pft_tag=""

  temp=$(mktemp -u --tmpdir=.)

  for npfts in $npfts_list
  do
    [[ $npfts -gt 0 ]] && pft_tag="${npfts}pfts_" || pft_tag=""
    outfile="${grid_name}_bc_land_frac_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    if [[ ${npfts} -gt 0 ]]
    then
      ${cdo} merge \
          ${grid_name}_notsea.nc \
          ${grid_name}_fract_glac.nc \
          ${grid_name}_sea.nc \
          ${grid_name}_fract_lake.nc \
          ${grid_name}_fract_land.nc \
          ${grid_name}_fract_veg.nc \
          ${grid_name}_veg_ratio_max_remap.nc \
          ${grid_name}_fract_pft*_remap.nc \
          ${grid_name}_land.nc \
          ${grid_name}_lake.nc \
          ${grid_name}_glac.nc \
          ${outfile}
    else
      ${cdo} merge \
          ${grid_name}_notsea.nc \
          ${grid_name}_fract_glac.nc \
          ${grid_name}_sea.nc \
          ${grid_name}_fract_lake.nc \
          ${grid_name}_fract_land.nc \
          ${grid_name}_fract_veg.nc \
          ${grid_name}_veg_ratio_max_remap.nc \
          ${grid_name}_land.nc \
          ${grid_name}_lake.nc \
          ${grid_name}_glac.nc \
          ${outfile}
    fi
    comment=$(cdo -s showattribute,comment ${input_file} | sed -e 's/\\n//g' -e 's/\"//g' | tr '\n' '; ' | tr -s " " | cut -d "=" -f2)
    ${cdo} setattribute,comment="${comment}" ${outfile} ${temp} && mv ${temp} ${outfile}

    #ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

  done
}

function modify_soil_and_root_depth {

  local soil_depth=${grid_name}_soil_depth_remap.nc
  local root_depth=${grid_name}_root_depth_remap.nc
  local maxmoist=${grid_name}_maxmoist_remap.nc
  local soil_field_cap=${grid_name}_soil_field_cap_remap.nc

  local temp=temp_$$

  echo "Modifying soil and root depth ..."

  # masks for tropics and extra-tropics (NH and SH)
  ${cdo} gtc,-1. -masklonlatbox,0.,360.,30.,90.   ${soil_depth} ${temp}_NH_mask.nc
  ${cdo} gtc,-1. -masklonlatbox,0.,360.,-90.,-30. ${soil_depth} ${temp}_SH_mask.nc
  ${cdo} gtc,-1. -masklonlatbox,0.,360.,-30.,30.  ${soil_depth} ${temp}_EQ_mask.nc

  # soil depth
  ${cdo} mulc,0.7 -mul ${soil_depth} -setmisstoc,0. ${temp}_EQ_mask.nc ${temp}_soil_depth_EQ.nc
  ${cdo} mul -addc,0.2 ${soil_depth} -setmisstoc,0. ${temp}_NH_mask.nc ${temp}_soil_depth_NH.nc
  ${cdo} mul -addc,0.2 ${soil_depth} -setmisstoc,0. ${temp}_SH_mask.nc ${temp}_soil_depth_SH.nc
  ${cdo} add ${temp}_soil_depth_EQ.nc ${temp}_soil_depth_NH.nc ${temp}_soil_depth_EQ_NH.nc
  ${cdo} add ${temp}_soil_depth_EQ_NH.nc ${temp}_soil_depth_SH.nc ${temp}_soil_depth_all.nc

  # compute root depth based on maximum root zone moisture and field capacity
  ${cdo} -setname,root_depth -div ${maxmoist} ${soil_field_cap} ${temp}_root_depth_all.nc

  # keep minimum values (0.3m for root depth and 0.5m for soil depth)
  ${cdo} gtc,0.5 ${temp}_soil_depth_all.nc ${temp}_soil_depth_all_05.nc
  ${cdo} setctomiss,0. ${temp}_soil_depth_all_05.nc ${temp}_soil_depth_all_05_miss.nc
  ${cdo} setmisstoc,0.5 -subc,1. ${temp}_soil_depth_all_05_miss.nc ${temp}_soil_depth_all_05_0.nc
  ${cdo} add ${temp}_soil_depth_all_05_0.nc -mul ${temp}_soil_depth_all.nc ${temp}_soil_depth_all_05.nc ${temp}_soil_depth_all_min.nc
  mv ${temp}_soil_depth_all_min.nc ${soil_depth}
  debug_file ${soil_depth}

  ${cdo} gtc,0.3 ${temp}_root_depth_all.nc ${temp}_root_depth_all_03.nc
  ${cdo} setctomiss,0. ${temp}_root_depth_all_03.nc ${temp}_root_depth_all_03_miss.nc
  ${cdo} setmisstoc,0.3 -subc,1. ${temp}_root_depth_all_03_miss.nc ${temp}_root_depth_all_03_0.nc
  ${cdo} add ${temp}_root_depth_all_03_0.nc -mul ${temp}_root_depth_all.nc ${temp}_root_depth_all_03.nc ${temp}_root_depth_all_min.nc
  mv ${temp}_root_depth_all_min.nc ${root_depth}
  debug_file ${root_depth}

  # clean
  rm ${temp}_*.nc

}

function debug_file {

  if [[ ${jsb4icon} == "true" ]]
  then
    local file=$1
    #local target_res=${source_res_dict[${bisection}]}
    local target_res=t255
    local debug_weights=wgt${debug_remap_scheme}_${grid_name}_to_${target_res}.nc

    if [[ ${debug} == "true" ]]
    then
      if [[ ! -f ${debug_weights} ]]
      then
        ${cdo} gen${debug_remap_scheme},${target_res}grid -random,${icon_grid_file} ${debug_weights}
      fi
      ${cdo} remap,${target_res}grid,${debug_weights} -setgrid,${icon_grid_file} -setmisstoc,0 ${file} debug_${file}
    fi
  fi
}
#_____________________________________________________________________________
#
# Main script
#_____________________________________________________________________________
#
mkdir -p ${workdir} || true
cd ${workdir}
echo "Workdir: ${workdir}"
rm -rf *

for bisection in $bisections $ggrids
do
  [[ $(echo ${bisection} | cut -c1) == T ]] && jsb4icon=false || jsb4icon=true
  if [[ ${jsb4icon} == true ]]; then
    bisect=$(echo ${bisection} | awk '{print $1 + 0}')
    if [[ $bisect -ge 8  ]]
    then
      # cdo="${cdo} -f nc4 -P 8"  # leads to crash with current cdo version on levante
      ulimit -s unlimited
      ulimit -v unlimited
      ulimit -m unlimited
    fi

    if [[ $bisect -le 5 ]]
    then
      remap_scheme=ycon
    else
      remap_scheme=dis
    fi

    grid_name=${refinement}
    source_res=${source_res_dict[${bisect}]}
  else         # jsbach4 on gaussian grid
    grid_name=${bisection}
    remap_scheme=ycon
    source_res=${grid_name}
    icon_grid=${grid_name}
  fi

  echo "=============================================================="
  echo "===  Generating files for $grid_name from $source_res      ==="
  echo "=============================================================="

  year=${start_year}

  if [[ ${jsb4icon} == true ]]; then
    icon_grid=icon_grid_${atmGridID}_${grid_name}_G.nc

    remap_weights=wgt${remap_scheme}_${source_res}_to_${grid_name}.nc


    case ${source_res} in
	  T* )
	    source_grid=${source_res}grid
	    ;;
	  05 )
	    ${cdo} griddes ${pool_prepare}/05/ECHAM6/05_jan_surf.nc > 0.5grid
	    source_grid=0.5grid
	    ;;
    esac

    command="${cdo} gen${remap_scheme},$icon_grid_file -setmissval,-9.e33 -random,${source_grid}  $remap_weights"
    echo "Generating remap weights: $command"
    time $command

    echo "extpar input file: $extpar_file"
  fi

  echo "==================================================================="
  echo "===  Get initial_tarfiles from git.mpimet.mpg.de/git/mpiesm.git ==="
  echo "==================================================================="

  if [[ ! -d initial_tarfiles ]]; then
    git clone -n -b mpiesm-landveg https://git.mpimet.mpg.de/git/mpiesm.git mpiesm-landveg
    # git clone -n -b mpiesm-landveg git@git.mpimet.mpg.de:mpiesm.git mpiesm-landveg
    cd mpiesm-landveg
    git config core.sparseCheckout true
    echo /contrib/initial_tarfiles >> .git/info/sparse-checkout
    git checkout f4a5f920a
    mv contrib/initial_tarfiles/ ..
    cd ..
    rm -rf mpiesm-landveg
  else
    echo " Initial_tarfiles already available"
  fi
  cd initial_tarfiles

  echo "==================================================================="
  echo Creating Gaussian input file ...
  generate_gauss_data

  cd ..

  shopt -s extglob
  #[[ $npfts -gt 0 ]] && \
  #    input_file=`ls ${input_root}/r0*/${source_res}/jsbach_${source_res}+(GR|TP)*_fractional_${npts}tiles_5layers_${year}.nc | tail -1` || \
  #    input_file=`ls ${input_root}/r0*/${source_res}/jsbach_${source_res}+(GR|TP)*_fractional_11tiles_5layers_${year}.nc | tail -1`
  input_file=`ls jsbach_${source_res}${res_oce}*_11tiles_5layers_${year}_no-dynveg.nc`
  shopt -u extglob

  echo "Extracting Gaussian variables from ${input_file}"
  ${cdo} splitname ${input_file} gauss_${source_res}_
  if [[ ${echam_fractional} == "true" ]]; then
    ${cdo} setctomiss,0 gauss_${source_res}_slf.nc gauss_${source_res}_notsea_miss.nc
  else
    ${cdo} setctomiss,0 gauss_${source_res}_slm.nc gauss_${source_res}_notsea_miss.nc
  fi
  ${cdo} ifnotthen gauss_${source_res}_glac.nc gauss_${source_res}_notsea_miss.nc \
         gauss_${source_res}_non-glac-land.nc

  echo Processing $icon_grid ...

  if [[ ${jsb4icon} == true ]]
  then

    generate_land_sea_mask

    for var in $varlist
    do
      extrapolate_var $var
    done

    for var in $varlist
    do
      remap_var $var
    done

  else
    for var in $varlist
    do
      pseudo_remap_var $var
    done
  fi

  generate_fractions

  crosscheck_vars

  if [[ ${jsb4icon} == true ]]
  then
    modify_soil_and_root_depth
  fi

  generate_output_files

  if [[ ${dryrun} != "true" ]]
  then
    echo output directory ${output_dir}
    mkdir -p ${output_dir}
  fi

  if [[ ${dryrun} != "true" ]]
  then
    echo "Copying output files to ${output_dir}"
    [[ ${clean} == "true" ]] && cpmv="mv" || cpmv="cp -p"
    [[ ${jsb4icon} == true ]] && grid_tag="" || grid_tag=${grid_name}${res_oce}_
    for npfts in $npfts_list
    do
      [[ $npfts -gt 0 ]] && pft_tag="${npfts}pfts_" || pft_tag=""
      ${cpmv} ${grid_name}_bc_land_frac_${pft_tag}$year.nc ${output_dir}/bc_land_frac_${pft_tag}${grid_tag}$year.nc
    done
    ${cpmv} ${grid_name}_bc_land_phys_$year.nc ${output_dir}/bc_land_phys_${grid_tag}$year.nc
    ${cpmv} ${grid_name}_bc_land_soil_$year.nc ${output_dir}/bc_land_soil_${grid_tag}$year.nc
    ${cpmv} ${grid_name}_bc_land_sso_$year.nc ${output_dir}/bc_land_sso_${grid_tag}$year.nc
    ${cpmv} ${grid_name}_ic_land_soil_$year.nc ${output_dir}/ic_land_soil_${grid_tag}$year.nc

    if [[ ${coupled} == "true" ]] && [[ ${hd_file}x != 'x' ]]
    then
      cp -p ${hd_file} ${output_dir}/bc_land_hd.nc
    fi

  fi

  # -------------------------------------------------------
  # now we generate bc_land_frac files for additional years
  # -------------------------------------------------------

  for ((year=${start_year}+1;year<=${end_year};year++))
  do
    echo "now processing bc_land_fact for year " ${year}
    cd initial_tarfiles
    generate_gauss_data
    cd ..

    input_file=`ls jsbach_${source_res}${res_oce}*_11tiles_5layers_${year}_no-dynveg.nc`
    echo "Extracting Gaussian variables from ${input_file}"

    ${cdo} splitname ${input_file} gauss_${source_res}_
    if [[ ${echam_fractional} == "true" ]]; then
      ${cdo} setctomiss,0 gauss_${source_res}_slf.nc gauss_${source_res}_notsea_miss.nc
    else
      ${cdo} setctomiss,0 gauss_${source_res}_slm.nc gauss_${source_res}_notsea_miss.nc
    fi

    if [[ ${jsb4icon} == true ]]
    then
      extrapolate_var   cover_fract
      remap_var         cover_fract
    else
      pseudo_remap_var  cover_fract
    fi

    generate_fractions
    crosscheck_vars

    generate_frac_output_files

    if [[ ${dryrun} != "true" ]]
    then
      echo "Copying output files to ${output_dir}"
      [[ ${clean} == "true" ]] && cpmv="mv" || cpmv="cp -p"
      for npfts in $npfts_list
      do
        if [[ $npfts -gt 0 ]]
        then
          pft_tag="${npfts}pfts_"
          ${cpmv} ${grid_name}_bc_land_frac_${pft_tag}$year.nc ${output_dir}/bc_land_frac_${pft_tag}$year.nc
        fi
      done
    fi
  done

done

if [[ ${clean} == "true" && ${workdir} != ${output_dir} ]]
then
  echo "Cleaning up ..."
  cd ${cwddir}
  rm -rf ${workdir}
fi

echo "==================================================================="
echo " $0 completed"
echo "==================================================================="

