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

#############################################################################################################################
### This Programm converts the landcover map files (LUH2) readable by JSBACH3 to corresponding files readable by JSBACH4. ###
### The Program  gauss_to_icon.ksh  in  jsbach/scripts/preprocessing/attic/  provides other infiles for the model.        ###
### The Program  jsbach4_ini_files_from_gauss.sh  in  jsbach/scripts/preprocessing/  generates from existing ECHAM6       ###
### files initial fields and boundary conditions for the land model JSBACH4.                                              ###
#############################################################################################################################

###############################################
# Some informations about what are we doing ###
###############################################

#   In the LUH2-file "cover_fract_T63_11tiles_1992.nc" the following structure of the grid box is assumed:
#                 box
#                 /  \
#          glacier    (vegetation)
#                     /    |    \
#                    pft1 pft2  pft3
#
#   In JS4 the following tile structure of the grid box is assumed:
#                 box
#                 /  \
#             lake    land
#                    /    \
#               glacier    vegetation
#                          /    |    \
#                        pft1 pft2  pft3
#
# Note, in the LUH2 file:
# * The cover_fracts sum up to 1 on the grid box. Not only for the 11 pfts but also over glacier, ocean or see areas!
# * The cover_fracts are relative to the grid box - as there are no parent tiles.
# * Glacier and pfts are disjunct. Thus a gridbox has either a glacier fract of 1 or 0.
# * There are no lake fractions. Instead there are pfts (even on the Caspian and Black Sea).
# * At the borders of the continents we have missing values (9.9e+36)!
# * Glacier cover_fracts are together with tropical trees on the 1. tile! The differentiation is given by the 1. cover_type.
#   C3 crops und C4 crops cover_fracts are together on the 11. tile! The differentiation is given by the 11. cover_type.
# * The cover_types give the lct numbers. In JS4 we do not need cover_types anymore, because we do not want two tiles on one record and
#   because the assignment to a lct number is done in the "usecase" and not in an input file.
#
# The LUH2 file includes the following pfts:
# Level  1: glacier and Tropical evergreen trees
# Level  2: Tropical deciduous trees
# Level  3: Extra-tropical evergreen trees
# Level  4: Extra-tropical deciduous trees
# Level  5: Raingreen shrubs
# Level  6: Deciduous shrubs
# Level  7: C3 grass
# Level  8: C4 grass
# Level  9: C3 pasture
# Level 10: C4 pasture
# Level 11: C3 and C4 crop

######################
### User Interface ###
######################
CDO="cdo -b F64"    # Set precision of cdo, e.g. "cdo -b F64"
OUT_DIR="Map-files_from-JS3-to-JS4"
INFILES_GRID="T63"  # e.g. T63
INFILES_PATH="/pool/data/JSBACH/input/r0010/${INFILES_GRID}/land_use_maps/" # The LUH2 data can be found here: /pool/data/JSBACH/input/r0010
                                                                            # E.g. /pool/data/JSBACH/input/r0010/${INFILES_GRID}/land_use_maps/
AUTO_LIST="yes" # Generate list of LUH2-MAPs-infiles automatically instead of naming the files explicitly.
# Y_START="1850"  # Start year of LUH2-MAPs-files that will be converted. cover_fract_T63_11tiles_${YY}.nc. Only used if AUTO_LIST="yes"
Y_START="1979"  # Start year of LUH2-MAPs-files that will be converted. cover_fract_T63_11tiles_${YY}.nc. Only used if AUTO_LIST="yes"
Y_END="2014"    # End year of LUH2-MAPs-files that will be converted. cover_fract_T63_11tiles_${YY}.nc. Only used if AUTO_LIST="yes"
#LIST_INFILE_NAME="cover_fract_T63_11tiles_1850.nc"   # Only used if AUTO_LIST="no". "Map"-infiles containing the cover fractions. E.g. "cover_fract_T63_11tiles_1976.nc"



MAKE_JS4_FRACTIONS_FILE="yes" # yes = the input map file are converted in a JS4 useable format. This does not mean that their grid is changed!
                              #       Only the cover_fracts are converted. This means C3 and C4 crops and tropical trees and glacier tiles are
                              #       splitted and put in separate records. Furthermore the variable names are adapted that JS4 can read them.
MAPS_LSM_PATH="/pool/data/JSBACH/input/r0010/T63" # Land-Sea-Mask that is used by the infiles.
                                                  # E.g. /pool/data/JSBACH/input/r0010/T63/
MAPS_LSM_FILE="bc_land_frac_11pfts_${INFILES_GRID}GR15_1976.nc" # The lsm-file that was used for the LUH2 maps. This is not necessarily the same
                                                                # as in the BC_LAND_FRAC_SOURCE_FILE.
                                                                # E.g. bc_land_frac_11pfts_T63GR15_1976.nc
SEPARATE_CROPS="no"    # yes for separating C3 and C4 crops as given in cover_type of the input file. Standard is no.
NUMBER_OF_PFTS="11"    # E.g. 11.



MAKE_REMAPPED_FILE="yes" # yes = convert maps to the grid given in NEW_GRID_FILE. E.g. R02B04.
BC_LAND_FRAC_SOURCE_FILE="/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/preliminary_land/r0002/R02B04_G/land/bc_land_frac_11pfts_1976.nc"
                         # File that defines the LSM of the remapped output files including possible additional grid cells after remapping the maps-files
                         # Nevertheless, this file should use the grid defined under NEW_GRID_FILE below
                         # coupled ocean:   e.g. /pool/data/ICON/coupled_input_temp/aloy_pre04.ff/Land/r0002/bc_land_frac_11pfts_1976.nc
                         # uncoupled ocean: e.g. /pool/data/ICON/grids/private/mpim/icon_preprocessing/source/preliminary_land/r0002/R02B04_G/land/bc_land_frac_11pfts_1976.nc
T255_DEBUG_FILES="no"    # yes = t255grid files for debugging are created. You can use them to regard the fields with ncview. Time-consuming.
                         # Note, differences may also occur from the remapping!
SOURCE_GRID="${INFILES_GRID}grid"   # For remapping: Name of the Grid of the LUH2 Map-Infiles
NEW_GRID_FILE="/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/grids/icon_grid_0005_R02B04_G.nc"                       # for AMIP
# NEW_GRID_FILE="/pool/data/ICON/coupled_input_temp/aloy_pre04.ff/ATMOOCEANINP_pre04_OceWithCoast_158km_editSLOHH2017_G.cdo172.nc"    # for old coupled runs aloi_pre04.ff
# NEW_GRID_FILE="/pool/data/ICON/grids/private/rene/mpim/0013/icon_grid_0013_R02B04_G.nc"                                               # for new coupled ruby runs
            # For remapping: Grid file that defines the grid for the remapped maps-files
            # E.g. /pool/data/ICON/coupled_input_temp/aloy_pre04.ff/ATMOOCEANINP_pre04_OceWithCoast_158km_editSLOHH2017_G.cdo172.nc  for R2B4 for ocean coupled run
            # E.g. /pool/data/ICON/grids/private/mpim/icon_preprocessing/source/grids/icon_grid_0005_R02B04_G.nc                     for R2B4 for ocean uncoupled run
OUT_GRID_NAME="R2B4" # Just for the filenames after remapping. E.g. R2B4
REMAPPING="SCRIP"    # NN for nearest neighbor else SRIP remapping is used. Results are very similar. However, SCRIP remapping is recommended.



MAKE_BC_LAND_FRAC_FILE="yes" # yes = create a JS4 readable and consistent bc_land_fract-file.
                            #       This means include the tiles notsea, sea, fract_lake, fract_land, veg_ratio_max, land, lake, glac, veg
BC_LAND_FRAC_OUT_FILE="bc_land_frac_11pfts" # This is only the base name, the base name of the input file will be added.



CHECK="no"  # yes = create a quick check file

##########################
### Make output folder ###
##########################

if [ ! -d ${OUT_DIR} ] ; then
   mkdir ${OUT_DIR}
   chmod 700 ${OUT_DIR}
else
   echo "There is already an output folder ${OUT_DIR}."
   echo "Overwrite it?"
   echo "YES =overwrite it, delete existing folder"
   echo "zzz =every other input keeps the folder and writes in it"
   read y
   if [ "${y}" == "YES" ] ; then
      echo "YES used, therefore delete folder..."
      rm -r ${OUT_DIR}
      mkdir ${OUT_DIR}
      chmod 700 ${OUT_DIR}
   fi
fi
HIER=$(pwd)
cd ${OUT_DIR}
mkdir -p tempfiles

################################
### Interpret interface data ###
################################

if [ ${AUTO_LIST} == "yes" ] ; then
  for YY in $(seq ${Y_START} ${Y_END});   do   # if a infile list shall be generated instead of explicitly name them
    YY=$(printf "%04d" ${YY})
    LIST_INFILE_NAME="${LIST_INFILE_NAME} cover_fract_${INFILES_GRID}_11tiles_${YY}.nc"
  done
fi

##########
### Go ###
##########

for INFILE_NAME in ${LIST_INFILE_NAME} ; do

  INFILE_NAME_BASE=$(echo ${INFILE_NAME} | cut -d. -f1 )

  # Convert the JS3 format "convention" to JS4 format
  if [ ${MAKE_JS4_FRACTIONS_FILE} == "yes" ] ; then
     # Prepare infile
     cp ${INFILES_PATH}/${INFILE_NAME}         .
     ${CDO}  -setmissval,-9.e33        ${INFILE_NAME}         ZWERG1
     ${CDO}  -delvar,soil_layer_depth  ZWERG1                 ZWERG2

     # Get rid of the fractions on the oceans
     cp ${MAPS_LSM_PATH}/${MAPS_LSM_FILE}      .
     ${CDO}  -selvar,notsea            ${MAPS_LSM_FILE}       MAPS_NOTSEA.nc
     ${CDO}  -mul                      ZWERG2                 MAPS_NOTSEA.nc     INFILE_MAPS.nc

     # Split into level=tile files
     ${CDO} -splitlevel  INFILE_MAPS.nc   MAPS_L_

     # Replace differentiation of tiles (levels to variable names) and split into single records
     for LL in $(seq 1 ${NUMBER_OF_PFTS});  do
        LL=$(printf "%02d" ${LL})
        echo ${LL}
        ${CDO}  -chname,cover_fract,fract_pft${LL}  -chname,cover_type,cover_type_pft${LL}   -setlevel,0     MAPS_L_0000${LL}.nc       MAPS_L_renamed_0000${LL}.nc
        ${CDO}  -splitrec  MAPS_L_renamed_0000${LL}.nc          MAPS_L_renamed_0000${LL}_
        mv                 MAPS_L_renamed_0000${LL}_000001.nc   MAPS_L_renamed_0000${LL}_cover_fract
        mv                 MAPS_L_renamed_0000${LL}_000002.nc   MAPS_L_renamed_0000${LL}_cover_type
     done

     # Split glacier and pft01
     ${CDO}  -eqc,1      MAPS_L_renamed_000001_cover_type  MAPS_L_renamed_000001_cover_type_glacier
     ${CDO}  -eqc,2      MAPS_L_renamed_000001_cover_type  MAPS_L_renamed_000001_cover_type_pft01
     ${CDO}  -chname,fract_pft01,fract_glac  -setmisstoc,0.000000   -ifthen  MAPS_L_renamed_000001_cover_type_glacier   MAPS_L_renamed_000001_cover_fract    MAPS_GLACIER_FRACT.nc
     ${CDO}  -setmisstoc,0.000000   -ifthen  MAPS_L_renamed_000001_cover_type_pft01     MAPS_L_renamed_000001_cover_fract    MAPS_PFT01_FRACT.nc

     # Split C3 and C4 crops
     if [ "${SEPARATE_CROPS}" == "yes" ] ; then
       ${CDO}  -eqc,20      MAPS_L_renamed_000011_cover_type  MAPS_L_renamed_000011_cover_type_CROP-C3
       ${CDO}  -eqc,21      MAPS_L_renamed_000011_cover_type  MAPS_L_renamed_000011_cover_type_CROP-C4
       ${CDO}  -chname,fract_pft11,fract_pft12 -setmisstoc,0.000000   -ifthen  MAPS_L_renamed_000011_cover_type_CROP-C4   MAPS_L_renamed_000011_cover_fract    CROP-C4_FRACT.nc
       ${CDO}                                              -setmisstoc,0.000000   -ifthen  MAPS_L_renamed_000011_cover_type_CROP-C3   MAPS_L_renamed_000011_cover_fract    CROP-C3_FRACT.nc
     fi

     # Create new fraction file in a JS4 readable format
     if [ "${SEPARATE_CROPS}" == "yes" ] ; then
        ${CDO}  -merge  MAPS_PFT01_FRACT.nc    MAPS_L_renamed_000002_cover_fract  MAPS_L_renamed_000003_cover_fract  MAPS_L_renamed_000004_cover_fract  MAPS_L_renamed_000005_cover_fract  MAPS_L_renamed_000006_cover_fract  MAPS_L_renamed_000007_cover_fract  MAPS_L_renamed_000008_cover_fract  MAPS_L_renamed_000009_cover_fract  MAPS_L_renamed_000010_cover_fract    CROP-C3_FRACT.nc  CROP-C4_FRACT.nc  MAPS_GLACIER_FRACT.nc        MAPS_${INFILE_NAME_BASE}_JS4-format.nc
     else
        ${CDO}  -merge  MAPS_PFT01_FRACT.nc  MAPS_L_renamed_000002_cover_fract  MAPS_L_renamed_000003_cover_fract  MAPS_L_renamed_000004_cover_fract  MAPS_L_renamed_000005_cover_fract  MAPS_L_renamed_000006_cover_fract  MAPS_L_renamed_000007_cover_fract  MAPS_L_renamed_000008_cover_fract  MAPS_L_renamed_000009_cover_fract  MAPS_L_renamed_000010_cover_fract  MAPS_L_renamed_000011_cover_fract   MAPS_GLACIER_FRACT.nc     MAPS_${INFILE_NAME_BASE}_JS4-format.nc
     fi

     # Sum all tiles up. Used later as land mask and also as a control file (should be 1 on each land grid box).
     ${CDO}  -expr,'summe=fract_pft01+fract_pft02+fract_pft03+fract_pft04+fract_pft05+fract_pft06+fract_pft07+fract_pft08+fract_pft09+fract_pft10+fract_pft11+fract_glac'    MAPS_${INFILE_NAME_BASE}_JS4-format.nc     MAPS_${INFILE_NAME_BASE}_JS4-FORMAT_SUM.nc

     # Clean up
     mv   INFILE_MAPS.nc  ZWERG1  ZWERG2  MAPS_L_*  MAPS_GLACIER_FRACT.nc   MAPS_PFT01_FRACT.nc   MAPS_NOTSEA.nc    tempfiles/

  fi


  # Remap the JS4-format-MAPS-file to the horizontal output grid of JS4.
  # The LSM of the LUH-MAPS-files do (normaly) not match with the LSM used by JS4.
  # As a consequence of the grids missmatch, grid cells may appear where we do not have values from the LUH-MAPS-files. There we put grass.
  if [ ${MAKE_REMAPPED_FILE} == "yes" ] ; then
    # Remap the JS4-format-MAP-file
    if [ ${REMAPPING} == "NN" ] ; then
      ${CDO}  -remapnn,${NEW_GRID_FILE}               MAPS_${INFILE_NAME_BASE}_JS4-format.nc        MAPS_${INFILE_NAME_BASE}_JS4-format_${OUT_GRID_NAME}.nc
    else
      ${CDO}  -genycon,${NEW_GRID_FILE}    -setmissval,-9.e33    -random,${SOURCE_GRID}   remap_weights
      ${CDO}  -remap,${NEW_GRID_FILE},remap_weights               MAPS_${INFILE_NAME_BASE}_JS4-format.nc        MAPS_${INFILE_NAME_BASE}_JS4-format_${OUT_GRID_NAME}.nc
    fi

    # Get land mask from JS4-format MAP-file
    if [ ${REMAPPING} == "NN" ] ; then
      ${CDO}  -remapnn,${NEW_GRID_FILE}               MAPS_${INFILE_NAME_BASE}_JS4-FORMAT_SUM.nc   ZWERG3
    else
      ${CDO}  -remap,${NEW_GRID_FILE},remap_weights   MAPS_${INFILE_NAME_BASE}_JS4-FORMAT_SUM.nc   ZWERG3
    fi
    # where pft-sum is not 1 set to zero.
    ${CDO}  -gec,0.999999999                          ZWERG3                                       MASK_${INFILE_NAME_BASE}_JS4-FORMAT_SUM_${OUT_GRID_NAME}.nc

    # Get land mask from bc_land_fract file (this LSM shall be used in the JS4 run)
    cp ${BC_LAND_FRAC_SOURCE_FILE}   .
    ${CDO}  -setmissval,-9.e33             ${BC_LAND_FRAC_SOURCE_FILE}  BC_INFILE.nc
    ${CDO}  -selvar,notsea                 BC_INFILE.nc                 BC_NOTSEA.nc
    ${CDO}  -gec,0.0001                    BC_NOTSEA.nc                 MASK_BC_NOTSEA.nc # costlines meight have fractions therefore this decision

    # Mask differences between the two land masks
    ${CDO}  -setctomiss,0.000000000   -sub   MASK_BC_NOTSEA.nc   MASK_${INFILE_NAME_BASE}_JS4-FORMAT_SUM_${OUT_GRID_NAME}.nc   MASK_BC_NOTSEA_minus_JS4-FORMAT_SUM_${OUT_GRID_NAME}.nc # here missing values exist in the file
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid    MASK_BC_NOTSEA_minus_JS4-FORMAT_SUM_${OUT_GRID_NAME}.nc     BC_NOTSEA_minus_JS4-FORMAT_SUM_T255.nc

    # Mask where ${BC_LAND_FRAC_SOURCE_FILE} notsea exists but in the JS4-format MAPS-file no land exists (there we will have to add values)
    ${CDO}  -setmisstoc,0.0000000000  -gec,0.000000001         MASK_BC_NOTSEA_minus_JS4-FORMAT_SUM_${OUT_GRID_NAME}.nc      MASK_ADDITIONAL_GP_${OUT_GRID_NAME}.nc # Note, there are missing values in the result
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid    MASK_ADDITIONAL_GP_${OUT_GRID_NAME}.nc                      MASK_ADDITIONAL_GP_T255.nc

    # Add the additional grid cells to the fract_pft08 (C3 grass)
    ${CDO}  -selvar,fract_pft08   MAPS_${INFILE_NAME_BASE}_JS4-format_${OUT_GRID_NAME}.nc   MAPS_fract_pft08_${OUT_GRID_NAME}.nc
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid    MAPS_fract_pft08_${OUT_GRID_NAME}.nc                  MAPS_fract_pft08_T255
    ${CDO}  -mul  MASK_BC_NOTSEA.nc                            MAPS_fract_pft08_${OUT_GRID_NAME}.nc    MAPS_fract_pft08_${OUT_GRID_NAME}.nc # -setctomiss,0.000000000 inserts miss. vals. for optical reasons
    ${CDO}  -max  MAPS_fract_pft08_${OUT_GRID_NAME}.nc   MASK_ADDITIONAL_GP_${OUT_GRID_NAME}.nc        NEW_fract_pft08_${OUT_GRID_NAME}.nc
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid   NEW_fract_pft08_${OUT_GRID_NAME}.nc                    NEW_fract_pft08_${OUT_GRID_NAME}_T255

    # Make sure that there are no values on the additional grid cells anywhere in the JS4-format MAP-file
    ${CDO}  -ifnotthen   MASK_ADDITIONAL_GP_${OUT_GRID_NAME}.nc   MAPS_${INFILE_NAME_BASE}_JS4-format_${OUT_GRID_NAME}.nc   MAPS_${INFILE_NAME_BASE}_JS4-format_NO_ADD_GP_${OUT_GRID_NAME}.nc

    # Now replace fract_pft08
    ${CDO}  -replace   MAPS_${INFILE_NAME_BASE}_JS4-format_NO_ADD_GP_${OUT_GRID_NAME}.nc    NEW_fract_pft08_${OUT_GRID_NAME}.nc  MAPS_${INFILE_NAME_BASE}_JS4-format_ADDITIONAL_GP_${OUT_GRID_NAME}.nc
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid    MAPS_${INFILE_NAME_BASE}_JS4-format_ADDITIONAL_GP_${OUT_GRID_NAME}.nc      MAPS_${INFILE_NAME_BASE}_JS4-format_ADDITIONAL_GP_T255.nc

    # Finally, we have to get rid of the dimension "ntiles" because JS4 can not read this
    module load nco
    ncks   -C -x -v ntiles      MAPS_${INFILE_NAME_BASE}_JS4-format_ADDITIONAL_GP_${OUT_GRID_NAME}.nc       ZWERG7
    ncwa   -O -a    ntiles      ZWERG7                                                                      ZWERG8
    rm ZWERG7

    # Sum all tiles up as a control file (should be 1 on each land grid box).
    ${CDO}  -expr,'summe=fract_pft01+fract_pft02+fract_pft03+fract_pft04+fract_pft05+fract_pft06+fract_pft07+fract_pft08+fract_pft09+fract_pft10+fract_pft11+fract_glac'    ZWERG8     MAPS_${INFILE_NAME_BASE}_JS4-format_NewSum_${OUT_GRID_NAME}.nc

    # Rename the remapped file
    if [ ${AUTO_LIST} == "yes" ] ; then
      YY=$(echo ${INFILE_NAME_BASE} | cut -d_ -f5 | cut -d. -f1)
      cp ZWERG8    REMAPPED_cover_fract_T63_11tiles_${YY}.nc
    else
      cp ZWERG8    REMAPPED_cover_fract_from_${INFILE_NAME_BASE}
    fi

   # Clean up
    mv  remap_weights  ZWERG3    MAPS_${INFILE_NAME_BASE}_JS4-FORMAT_SUM.nc   MAPS_${INFILE_NAME_BASE}_JS4-format_NewSum_${OUT_GRID_NAME}.nc   BC_NOTSEA.nc    tempfiles/
  fi


  # JS4 reads the cover fractions at the beginning from the bc_land_fract file. This file additionally includes the fractions for:
  # notsea, sea, fract_lake, fract_land, veg_ratio_max, land, lake, glac, fract_veg.
  # fract_veg has to be summed up from the new pft fractions. Note, in JS4 pft fractions are read in as relative to the veg tile.
  # The pft fractions sum um to 1 on the veg fraction. Glacier and veg fractions are disjunct. In JS4 fract_glac and fract_veg is used
  # realtive to the land and not to the box tile (as in the LUH2-infile)!
  # However, this is perfect. As we have lakes in JS4 we want to scale them down to the land tile.
  if [ ${MAKE_BC_LAND_FRAC_FILE} == "yes" ] ; then
    # Missing values should normaly be handeled by jsbach, but it is not in the moment
    ${CDO}  -setmisstoc,0.0000000000   ZWERG8    ZWERG9

    # Calculate new fract_veg
    ${CDO}  -expr,'fract_veg=fract_pft01+fract_pft02+fract_pft03+fract_pft04+fract_pft05+fract_pft06+fract_pft07+fract_pft08+fract_pft09+fract_pft10+fract_pft11'    ZWERG9    MAPS_${INFILE_NAME_BASE}_fract_veg_${OUT_GRID_NAME}.nc
              [ ${T255_DEBUG_FILES} == "yes" ] && ${CDO}  -remapycon,t255grid    MAPS_${INFILE_NAME_BASE}_fract_veg_${OUT_GRID_NAME}.nc     MAPS_${INFILE_NAME_BASE}_fract_veg_T255.nc

    # Now we can create our new bc_land_fract-file
    ${CDO}  -selvar,notsea,sea,fract_lake,fract_land,veg_ratio_max,land,lake,glac    BC_INFILE.nc   BC_INFILE_WITHOUT_VEG.nc
    # Attention: GLACIERS are disjunct with pfts (because of the MAPS-infile), but pfts and glaciers of the MAP-files do not necessarily agree with the BC_LAND_FRAC_SOURCE_FILE!
    ${CDO}  -merge   BC_INFILE_WITHOUT_VEG.nc   ZWERG9  MAPS_${INFILE_NAME_BASE}_fract_veg_${OUT_GRID_NAME}.nc   ZWERG4

    # Finally, we have to get rid of the dimension "ntiles" because JS4 can not read this
    # module load nco
    # ncks   -C -x -v ntiles      ZWERG4       ZWERG5
    # ncwa   -O -a    ntiles      ZWERG5       ${BC_LAND_FRAC_OUT_FILE}__from_${INFILE_NAME}
    # rm ZWERG5
    mv ZWERG4 ${BC_LAND_FRAC_OUT_FILE}__from_${INFILE_NAME}
    # Clean up
    # mv    ZWERG4     BC_INFILE_WITHOUT_VEG.nc   tempfiles/
  fi


  # Simple check by remapping pft 2 of the infile and/or the outfile. pft 2 is unchanged in the program, only remapped.
  # As the interpolation error should be larger at the borders of the pft the differences should be largest there
  if [ ${CHECK} == "yes" ] ; then
    # Remap both to t255grid
    ${CDO} -remapycon,t255grid -selvar,fract_pft02                                              ${BC_LAND_FRAC_OUT_FILE}__from_${INFILE_NAME}  CHECK_out_fract_pft02
    ${CDO} -remapcon,t255grid  -setrtoc,-0.0000001,0.0000001,0.0  -selvar,cover_fract   -sellevel,2   ${INFILE_NAME}                                 CHECK_in_cover_fract_level2
    ${CDO} -setrtomiss,-0.000001,0.000001  -sub    CHECK_out_fract_pft02    CHECK_in_cover_fract_level2   CHECK_diff_fract_pft02__from_${INFILE_NAME} # insert miss. val. to avoid missleading high percentage values in the next step
    ${CDO} -div  -mulc,100.0    CHECK_diff_fract_pft02__from_${INFILE_NAME}   CHECK_in_cover_fract_level2     CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc

#    # Remap infile to outgrid (interpolation error cancels out)
#    # As the interpolation error should be more or less the same for in- and outfile
#    ${CDO}                     -selvar,fract_pft02                                               ${BC_LAND_FRAC_OUT_FILE}__from_${INFILE_NAME}   CHECK_out_fract_pft02
#    if [ ${REMAPPING} == "NN" ] ; then
#      ${CDO} -remapnn,${NEW_GRID_FILE}                         -setrtoc,-0.0000001,0.0000001,0.0  -selvar,cover_fract   -sellevel,2   ${INFILE_NAME}      CHECK_in_cover_fract_level2
#    else
#      #${CDO} -remap,${NEW_GRID_FILE},tempfiles/remap_weights  -setrtoc,-0.0000001,0.0000001,0.0  -selvar,cover_fract   -sellevel,2   ${INFILE_NAME}      CHECK_in_cover_fract_level2
#     fi
#    ${CDO} -setrtomiss,-0.000001,0.000001  -sub    CHECK_out_fract_pft02    CHECK_in_cover_fract_level2   CHECK_diff_fract_pft02__from_${INFILE_NAME} # insert miss. val. to avoid missleading high percentage values in the next step
#    ${CDO} -div  -mulc,100.0    CHECK_diff_fract_pft02__from_${INFILE_NAME}   CHECK_in_cover_fract_level2     CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc
#    #${CDO} -infov  CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc
#    #${CDO} -remapycon,t255grid  CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc  CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc_t255

    # Remap the outfile back to the input grid e.g. to t63grid
#    ${CDO} -remapycon,${INFILES_GRID}grid -selvar,fract_pft02              ${BC_LAND_FRAC_OUT_FILE}__from_${INFILE_NAME}  CHECK_out_fract_pft02
#    ${CDO} -setrtoc,-0.0000001,0.0000001,0.0  -selvar,cover_fract   -sellevel,2  ${INFILE_NAME}                                 CHECK_in_cover_fract_level2
#    ${CDO} -sub    CHECK_out_fract_pft02    CHECK_in_cover_fract_level2   CHECK_diff_fract_pft02__from_${INFILE_NAME}
#    ${CDO} -div  -mulc,100.0    CHECK_diff_fract_pft02__from_${INFILE_NAME}   CHECK_in_cover_fract_level2     CHECK_diff_fract_pft02__from_${INFILE_NAME}_in_percent.nc


    # Clean up
    mv  CHECK_out_fract_pft02  CHECK_in_cover_fract_level2    tempfiles/
  fi
done

######################
### Finish message ###
######################

echo -e "\033[01;31;40m"
echo "Reached normal end of script!"
tput sgr0
exit
