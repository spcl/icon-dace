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
#SBATCH --job-name=extpar
#SBATCH --partition=compute
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --account=mj0060
#SBATCH --output=extpar.o%j
#SBATCH --error=extpar.e%j
#_____________________________________________________________________________
set -e
ulimit -s 262144
ulimit -c 0

echo "------------------------------------------------------------------------"
echo "-- Running $0"
echo "------------------------------------------------------------------------"

icon_grid=${icon_grid:-/pool/data/ICON/grids/public/mpim/0043/icon_grid_0043_R02B04_G.nc}

tile_mode=0
urban_mode=0

bindir=${extpar_dir:-/work/mj0060/m212070/test_extpar4jsbach/extpar4jsbach}/bin
pydir=${extpar_dir:-/work/mj0060/m212070/test_extpar4jsbach/extpar4jsbach}/python

workdir=${work_dir:-/work/mj0060/m212070/test_extpar4jsbach/extpar4jsbach/work}
extpar_file=${extpar_file:-icon_extpar__20230404.nc}

datadir=/work/mj0060/m212005/data/extpar-input-data

clean_up=true
cdo="cdo -s"

if [[ ! -d ${workdir} ]] ; then
  mkdir -p ${workdir}
fi
cd ${workdir}
echo "$0: Working directory: $(pwd)"

griddir=${icon_grid%/*}
icon_grid_file=${icon_grid##*/}  # i.e. ${atmGridID}_${refinement}_G

export OMP_NUM_THREADS=1
export NETCDF_OUTPUT_FILETYPE=NETCDF4

export PYTHONPATH=$PYTHONPATH:${pydir}/lib

if [ "${tile_mode}" = 1 ]
then
    echo "$0: TILE_MODE for EXTPAR is switched ON"
else
    echo "$0: TILE_MODE for EXTPAR is switched OFF"
fi
if [ "${urban_mode}" = 1 ]
then
    echo "$0: URBAN_MODE for EXTPAR is switched ON"
else
    echo "$0: URBAN_MODE for EXTPAR is switched OFF"
fi

#nest_TAG=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 1)
#nest_ID=$(basename ${icon_grid_file} .nc | awk -F '_' '{print $5}' | cut -c 3)
#echo "Check to recognize nested domain with nest_TAG=${nest_TAG} and nest_ID=${nest_ID}"

binary_extpar_consistency_check=extpar_consistency_check.exe
binary_aot=extpar_aot_to_buffer.exe
binary_cru=extpar_cru_to_buffer.py
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_ndvi=extpar_ndvi_to_buffer.py
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_alb=extpar_alb_to_buffer.py
binary_era=extpar_era_to_buffer.py

cat > INPUT_grid_org << EOF
&GRID_DEF
 igrid_type = 1,
 domain_def_namelist='INPUT_ICON_GRID'
/
EOF

cat > INPUT_ICON_GRID << EOF
&icon_grid_info
  icon_grid_dir='${griddir}'
  icon_grid_nc_file ='${icon_grid_file}'
/
EOF

if [ "${tile_mode}" = 1 ]
then
    grib_output_filename="${extpar_file%.nc}_tiles.g2"
    netcdf_output_filename="${extpar_file%.nc}_tiles.nc"
    echo "$0: Output file is $netcdf_output_filename"
else
    grib_output_filename="${extpar_file%.nc}.g2"
    netcdf_output_filename="${extpar_file}"
    echo "$0: Output file is $netcdf_output_filename"
fi

raw_data_alb='MODIS_alb.nc'
raw_data_alnid='MODIS_alnid.nc'
raw_data_aluvd='MODIS_aluvd.nc'
buffer_alb='albedo_buffer.nc'
output_alb='albedo_icon.nc'

raw_data_aot='aerosols/aot_GACP.nc'
buffer_aot='extpar_aot_BUFFER.nc'
output_aot='aot_extpar_ICON.nc'

raw_data_tclim_coarse='absolute_hadcrut3.nc'
raw_data_tclim_fine='CRU_T2M_SURF_clim.nc'
buffer_tclim='cru_buffer.nc'
output_tclim='crut_extpar_ICON.nc'

raw_data_glc2000='glc2000_byte.nc'
buffer_glc2000='extpar_landuse_BUFFER.nc'
output_glc2000='extpar_landuse_ICON.nc'
raw_data_glcc='landuse/GLCC_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_BUFFER.nc'
output_glcc='glcc_landuse_ICON.nc'

# raw_data_globcover='GLOBCOVER_L4_200901_200912_V2.3_int16.nc'
raw_data_globcover_0='/landuse/GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='/landuse/GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='/landuse/GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='/landuse/GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='/landuse/GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='/landuse/GLOBCOVER_5_16bit.nc'
buffer_lu='extpar_landuse_BUFFER.nc'
output_lu='extpar_landuse_ICON.nc'

raw_data_globe_A10='/topo/globe/GLOBE_A10.nc'
raw_data_globe_B10='/topo/globe/GLOBE_B10.nc'
raw_data_globe_C10='/topo/globe/GLOBE_C10.nc'
raw_data_globe_D10='/topo/globe/GLOBE_D10.nc'
raw_data_globe_E10='/topo/globe/GLOBE_E10.nc'
raw_data_globe_F10='/topo/globe/GLOBE_F10.nc'
raw_data_globe_G10='/topo/globe/GLOBE_G10.nc'
raw_data_globe_H10='/topo/globe/GLOBE_H10.nc'
raw_data_globe_I10='/topo/globe/GLOBE_I10.nc'
raw_data_globe_J10='/topo/globe/GLOBE_J10.nc'
raw_data_globe_K10='/topo/globe/GLOBE_K10.nc'
raw_data_globe_L10='/topo/globe/GLOBE_L10.nc'
raw_data_globe_M10='/topo/globe/GLOBE_M10.nc'
raw_data_globe_N10='/topo/globe/GLOBE_N10.nc'
raw_data_globe_O10='/topo/globe/GLOBE_O10.nc'
raw_data_globe_P10='/topo/globe/GLOBE_P10.nc'

buffer_topo='topography_BUFFER.nc'
output_topo='topography_ICON.nc'

buffer_era='era_buffer.nc'

raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'
output_ndvi='ndvi_extpar_ICON.nc'

raw_data_soil_FAO='soil/FAO_DSMW_double.nc'
raw_data_soil_HWSD='soil/HWSD0_30_topsoil.nc'
raw_data_deep_soil='soil/HWSD30_100_subsoil.nc'
buffer_soil='SOIL_BUFFER.nc'
output_soil='SOIL_ICON.nc'

raw_lookup_table_HWSD='soil/LU_TAB_HWSD_UF.data'
raw_HWSD_data='soil/HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='soil/HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='soil/HWSD_DATA_COSMO_EXTPAR.asc'

raw_data_flake='flake/GLDB_lakedepth.nc'
buffer_flake='flake_BUFFER.nc'
output_flake='ext_par_flake_ICON.nc'

# NOAA imperviou surface area dataset (default)
raw_data_isa_0='NOAA_ISA_16bit.nc'

# # EEA impervious surface area dataset
# raw_data_isa_0='EEA_ISA_4_16bit.nc'


buffer_isa='ISA_BUFFER.nc'
output_isa='ISA_extpar_ICON.nc'

# this file is adapted from Flanner (2009)
raw_data_ahf='AHF_2006_2.5min_latreverse.nc'

# # AHF is redistributed at 25km scales according to NOAA ISA.
# raw_data_ahf='AHF_2006_NOAAISAredistr.nc'

buffer_ahf='AHF_BUFFER.nc'
output_ahf='AHF_extpar_ICON.nc'

# create input namelists
cat > INPUT_AOT << EOF
&aerosol_raw_data
  raw_data_aot_path='${datadir}',
  raw_data_aot_filename='${raw_data_aot}'
/
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}'
/
EOF
#---
cat >> ${pydir}/lib/namelist.py << EOF
input_tclim = {
    'raw_data_t_clim_path': '${datadir}/cru/',
    'raw_data_tclim_coarse': '${raw_data_tclim_coarse}',
    'raw_data_tclim_fine': '${raw_data_tclim_fine}',
    't_clim_buffer_file': 'crutemp_climC_extpar_BUFFER.nc',
    'it_cl_type': 1
}
EOF

#---
cat > INPUT_LU << EOF
&lu_raw_data
   raw_data_lu_path='',
   raw_data_lu_filename='${datadir}${raw_data_globcover_0}' '${datadir}${raw_data_globcover_1}' '${datadir}${raw_data_globcover_2}' '${datadir}${raw_data_globcover_3}' '${datadir}${raw_data_globcover_4}' '${datadir}${raw_data_globcover_5}',
   i_landuse_data=1,
   ilookup_table_lu=1,
   ntiles_globcover=6
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}'
/
&glcc_raw_data
   raw_data_glcc_path='${datadir}',
   raw_data_glcc_filename='${raw_data_glcc}'
/
&glcc_io_extpar
   glcc_buffer_file='${buffer_glcc}'
/
EOF
#---
cat > INPUT_ORO << EOF
&oro_runcontrol
  lcompute_sgsl=.FALSE.
/
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type=1,
 lsso_param=.TRUE.,
 lsubtract_mean_slope=.FALSE.,
 raw_data_orography_path='',
 ntiles_column=4,
 ntiles_row=4,
 topo_files='${datadir}${raw_data_globe_A10}' '${datadir}${raw_data_globe_B10}' '${datadir}${raw_data_globe_C10}' '${datadir}${raw_data_globe_D10}' '${datadir}${raw_data_globe_E10}' '${datadir}${raw_data_globe_F10}' '${datadir}${raw_data_globe_G10}' '${datadir}${raw_data_globe_H10}' '${datadir}${raw_data_globe_I10}' '${datadir}${raw_data_globe_J10}' '${datadir}${raw_data_globe_K10}' '${datadir}${raw_data_globe_L10}' '${datadir}${raw_data_globe_M10}' '${datadir}${raw_data_globe_N10}' '${datadir}${raw_data_globe_O10}' '${datadir}${raw_data_globe_P10}'
/
EOF
cat > INPUT_OROSMOOTH << EOF
&orography_smoothing
  lfilter_oro=.FALSE.   ! Smoothing does not work for ICON grids
/
EOF
#---
cat > INPUT_RADTOPO << EOF
&radtopo
  lradtopo=.FALSE.,
  nhori=24,
/
EOF
#---
cat >> ${pydir}/lib/namelist.py << EOF
input_era = {
    'iera_type': 1,
    'raw_data_era_path': '${datadir}/era/',
    'raw_data_era_ORO': 'ERA5_ORO_1990.nc',
    'raw_data_era_T2M': 'ERA5_T2M_1990_2019.nc',
    'raw_data_era_SST': 'ERA5_SST_1990_2019.nc',
    'raw_data_era_SD': 'ERA5_SD_1990_2019.nc',
    'era_buffer_file': '${buffer_era}',
}
EOF
#---
cat >> ${pydir}/lib/namelist.py << EOF
input_ndvi = {
    'raw_data_ndvi_path': '${datadir}/ndvi',
    'raw_data_ndvi_filename': '${raw_data_ndvi}',
    'ndvi_buffer_file': '${buffer_ndvi}',
    'ndvi_output_file': '${output_ndvi}'
}
EOF
#---
cat > INPUT_SOIL << EOF
&soil_raw_data
 isoil_data = 2,
 ldeep_soil = .true.,
 raw_data_soil_path='${datadir}',
 raw_data_soil_filename='${raw_data_soil_HWSD}'
 raw_data_deep_soil_filename='${raw_data_deep_soil}'
/
&soil_io_extpar
  soil_buffer_file='${buffer_soil}',
  soil_output_file_consistent='${output_soil}'
/
&HWSD_index_files
 path_HWSD_index_files='${datadir}',
 lookup_table_HWSD='${raw_lookup_table_HWSD}',
 HWSD_data='${raw_HWSD_data}',
 HWSD_data_deep='${raw_HWSD_data_deep}',
/
EOF
#---
if [ "${urban_mode}" = 1 ]
then
    cat > INPUT_ISA << EOF
&isa_raw_data
   raw_data_isa_path='',
   raw_data_isa_filename='${datadir}${raw_data_isa_0}'
   ntiles_isa=1
/
&isa_io_extpar
   isa_buffer_file='${buffer_isa}',
   isa_output_file='${output_isa}'
/
EOF
    #---
    cat > INPUT_AHF << EOF
&ahf_raw_data
  raw_data_ahf_path='${datadir}',
  raw_data_ahf_filename='${raw_data_ahf}'
/
&ahf_io_extpar
 ahf_buffer_file='${buffer_ahf}',
 ahf_output_file='${output_ahf}'
/
EOF
fi
#---
cat > INPUT_FLAKE << EOF
&flake_raw_data
   raw_data_flake_path='${datadir}',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
/
EOF
#   flake_output_file='${output_flake}'
#---
cat >> ${pydir}/lib/namelist.py << EOF
input_alb = {
    'ialb_type': 1,
    'raw_data_alb_path': '${datadir}/albedo/',
    'raw_data_alb_filename': '${raw_data_alb}',
    'raw_data_alnid_filename': '${raw_data_alnid}',
    'raw_data_aluvd_filename': '${raw_data_aluvd}',
    'alb_buffer_file': '${buffer_alb}',
    'alb_output_file': '${output_alb}',
}
EOF
#_______________________________________________________________________________
# consistency check
cat > INPUT_CHECK << EOF
&extpar_consistency_check_io
  grib_output_filename='${grib_output_filename}'
  netcdf_output_filename='${netcdf_output_filename}'
  i_lsm_data=1
  land_sea_mask_file=""
  number_special_points=0
  tile_mode=${tile_mode}
/
EOF

#---
cat > INPUT_TCLIM << EOF
&t_clim_raw_data
  raw_data_t_clim_path='${datadir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
/

&t_clim_io_extpar
  t_clim_buffer_file='${buffer_tclim}',
  t_clim_output_file='${output_tclim}'
/
EOF
if [ "${urban_mode}" = 1 ]
then
    ${bindir}/${binary_ahf}
    ${bindir}/${binary_isa}
fi

set -x
${bindir}/${binary_topo}
${bindir}/${binary_soil}
${bindir}/${binary_alb}
${bindir}/${binary_ndvi}
${bindir}/${binary_flake}
${bindir}/${binary_lu}
${bindir}/${binary_cru}
${bindir}/${binary_aot}
${bindir}/${binary_era}
set +x

# the consistency check requires the output of
#   ${binary_aot}, ${binary_cru}, ${binary_lu}, ${binary_globe},
#   ${binary_ndvi}, ${binary_soil} and ${binary_flake}

cat > INPUT_TCLIM_FINAL << EOF
&t_clim_raw_data
  raw_data_t_clim_path='${datadir}',
  raw_data_t_clim_filename='${raw_data_tclim_fine}',
/

&t_clim_io_extpar
  t_clim_buffer_file='crutemp_climF_extpar_BUFFER.nc',
  t_clim_output_file='crutemp_climC_extpar_BUFFER.nc'
/
EOF
${bindir}/${binary_extpar_consistency_check}

# Final corrections for the Extpar output file
#   - Use SOILTYPE_2 instead of SOILTYP (zero everywhere)
#   - Set the ICON grid (clon and clat first removed to avoid conflict)
${cdo} delvar,SOILTYP,clon,clat $workdir/${netcdf_output_filename} $workdir/${netcdf_output_filename}.tmp
${cdo} chname,SOILTYP_2,SOILTYP $workdir/${netcdf_output_filename}.tmp $workdir/${netcdf_output_filename}.tmp2

#    - Use the glacier fraction from ${binary_lu} output, as the glacier mask correction
#      within the consistency check is not helpful in coupled setups, where the model land
#      sea mask does not match FR_LAND.
${cdo} setvar,ICE     -setlevel,0 -selvar,ICE_GLCC     glcc_landuse_BUFFER.nc ICE_GLCC.nc
${cdo} setvar,FR_LAND -setlevel,0 -selvar,FR_LAND_GLCC glcc_landuse_BUFFER.nc FR_LAND_GLCC.nc
${cdo} -O merge ICE_GLCC.nc FR_LAND_GLCC.nc  glcc_lu_vars.nc
${cdo} replace ${netcdf_output_filename}.tmp2 glcc_lu_vars.nc ${netcdf_output_filename}.tmp3

${cdo} setgrid,${griddir}/${icon_grid_file} $workdir/${netcdf_output_filename}.tmp3 \
    ${extpar_dir}/${netcdf_output_filename}

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "External parameters for ICON-Land model generated:"
echo "      ${extpar_dir}/${netcdf_output_filename}"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# Clean up
if [[ ${clean_up} == true ]]; then
  cd ..
  rm -fr $workdir
fi

exit 0

#echo 'Generation of global NetCDF attributes' 
#
#sed -n '/Code information/,/code information/p' log_consistency_check.txt > attrib_comment.txt
#library_name=$(fgrep  "Library name" attrib_comment.txt | awk -F ":"  '{print $2}')
#tag_name=$(fgrep  "Tag name" attrib_comment.txt | awk -F ":"  '{print $2}')
#revision_number=$(fgrep  "Revision number" attrib_comment.txt | awk -F ":"  '{print $2}')
#checkin_date=$(fgrep  "Checkin-Date" attrib_comment.txt | awk -F ":"  '{print $2}')
#code_modified=$(fgrep  "Code is modified" attrib_comment.txt | awk -F ":"  '{print $2}')
#compile_date=$(fgrep  "Compile-Date" attrib_comment.txt | awk -F ":"  '{print $2}')
#compiled_by=$(fgrep  "Compiled by" attrib_comment.txt | awk -F ":"  '{print $2}')
#compiled_on=$(fgrep  "Compiled on" attrib_comment.txt | awk -F ":"  '{print $2}')
#start_time=$(fgrep  "Current start time" attrib_comment.txt | awk -F ":"  '{print $2}')
#grib_md5sum=$(md5sum -b ${grib_output_filename})
#
#echo "MD5SUM of the GRIB file is $grib_md5sum"
#
#if [ "${nest_TAG}" = "N" ]
#then
#    if [ "${nest_ID}" = "2" ]
#    then
#        /usr/local/pkg/grib_api/1.12.3/CRAY/bin/grib_set -sgeneratingProcessIdentifier=2 ${grib_output_filename} ${grib_output_filename}_GPI2
#        echo "For nested domain nest_TAG = ${nest_TAG} with generatingProcessIdentifier=${nest_ID} the result file in GRIB2 is ${grib_output_filename}_GPI2"
#    fi
#    if [ "${nest_ID}" = "3" ]
#    then
#        /usr/local/pkg/grib_api/1.12.3/CRAY/bin/grib_set -sgeneratingProcessIdentifier=3 ${grib_output_filename} ${grib_output_filename}_GPI3
#        echo "For nested domain nest_TAG = ${nest_TAG} with generatingProcessIdentifier=${nest_ID} the result file in GRIB2 is ${grib_output_filename}_GPI3"
#    fi
#fi
#
#/e/uhome/jhelmert/bin/nco-4.4.7/src/nco/ncatted -h -a library_name,global,o,c,"$library_name"\
# -a tag_name,global,o,c,"$tag_name"\
# -a comment,global,d,, \
# -a history,global,d,, \
# -a bin_dir,global,o,c,"$bindir"\
# -a revision_number,global,o,c,"$revision_number"\
# -a checkin_date,global,o,c,"$checkin_date"\
# -a code_modified,global,o,c,"$code_modified"\
# -a compile_date,global,o,c,"$compile_date"\
# -a compiled_by,global,o,c,"$compiled_by"\
# -a compiled_on,global,o,c,"$compiled_on"\
# -a start_time,global,o,c,"$start_time"\
# -a tile_mode,global,o,c,"$tile_mode"\
# -a grib_md5sum,global,o,c,"$grib_md5sum" ${netcdf_output_filename}
#
#/e/uhome/jhelmert/bin/nco-4.4.7/src/nco/ncbo -O $workdir/${netcdf_output_filename} ${testfile}  diff.nc; /e/uhome/jhelmert/bin/cdo_xce infov diff.nc
#
#echo "Test of external parameters against version ${testfile} has been performed"
#echo "Difference to former version using cdo_xce infov diff.nc should be min,max,avg=0 for every parameter!"

#echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#echo "External parameters for refinement grid of ICON model generated in $workdir"
#echo "$testfile has been used as reference. Please check $workdir/diff.nc"

#echo "Check consistency to former versions with ncdiff -O ${netcdf_output_filename} /e/uscratch/jhelmert/ICON_EXTPAR_20141202/icon_extpar_0026_R03B07_G_20141202.nc diff.nc; cdo infov diff.nc"

#echo "For plotting results with bplot use: export ICON_COORDINATE_FILE=$griddir/$icon_grid_file"
#echo "or modify ${HOME}/.profile_bplot  "
#echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

