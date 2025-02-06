#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
#==========================================================================
#      Driver script for postprocessing and visulization for the 
#==========================================================================
#
# History: 
# Initial version by Stephan Lorenz (MPI-M, 2011-01)
#
# Short description:
#  - first test version: 1 plot of elevation after 3 days = 12 plot_steps (6-h writing)
#  - update to exp.hom_lsm_flat for 1 day (2011-06)
#  - now using "planet" bathymetry, experiment name including bathy-file (2011-10)
#
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   module ncl/5.2.0-bin for visualization using icon_plot.ncl
#
#==============================================================================

check_error()
{

# Check if the first parameter (return status) is not OK 

  if [ $1 -ne 0 ] 
  then

# Stop running this script and return the error status
    echo "ERROR: $2"
    exit $1
  fi
}

set -ex
#
# Check if a parameter is given

if [ "x$1" != "x" ]
then
  set_env=$1 
else
  set_env=/null
fi

if [ -f ${set_env} ] 
then 
  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi

#==========================================================================
# The directory in which the model output can be found. 
# Don't forget the trailing "/".
#
# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
icon_path=${basedir}/
resolution=${horizontal_resolution}${vertical_resolution}
#
# Where should the plot files be located? Don't forget the trailing "/".
model_data_path="${icon_path}experiments/$EXP/"
ncl_script_path="${icon_path}scripts/postprocessing/tools"
#plot_file_path="${model_data_path}plots/"

#==========================================================================
echo
echo "**********************************************************"
echo "***       ICON Tool Kit for Idealized Test Cases       ***"
echo "**********************************************************"
echo 
echo "=== Postprocessing started for the hydrostatic ocean $EXP Shallow Water test."

mkdir -p $model_data_path
cd $model_data_path

# Parameters for ncl including double quotes
ifile="${ExpName}_0001.nc"
ofile="${ExpName}_${resolution}_ELEV"
ofile="${plotBaseName}_W"                       #  Vertical velocity
otype="eps"                                     #  plot file format
varname="W"                                     #  Vertical velocity
varname="ELEV"                                  #  Surface elevation
maxvar=8                                        #  initial elevation: +-10 m
minvar=-8
selmode="halflog"
selmode="manual"                                #  use max/minvar
selmode="auto"                                  #  no max/minvar used
MAP="'mapLLC=(/-60,-35/)' 'mapURC=(/60,35/)'"   #  use local lat/lon-map (at equator)
MAP="'mapLLC=(/-95,  5/)' 'mapURC=(/40,85/)'"   #  use local lat/lon-map (10-80N)
MAP=""                                          #  global map
plotstep=12                                     #  3 days, 6-hourly data
plotstep=4                                      #  1 day,  6-hourly data
mname="wet_c"                                   #  grey shading outside basin
PROJSAT="projSat=True"                          #  Projection: Sattelite view
PROJSAT=" "                                     #  Projection: Cylindrical Equidistant
bstrg="S.J.Lorenz;"                             #  base string
tstrg="ICON Ocean"                              #  title string

# set the CDO path for the dwd to a local installation
case "$(hostname -d)" in
  dwd.de)
    CDO="/e/uhome/extrmuel/local/bin/cdo"
    ;;
  *)
    CDO="cdo"
esac
# these parameters are stored in here-document due to necessary single and double quotes:
cat >scr_${EXP}_nclcmd.here <<eo_here
  ncl ${ncl_script_path}/icon_plot.ncl \
      'iFile="$ifile"' 'oFile="$ofile"' 'varName="$varname"' 'oType="$otype"' 'altLibDir="${ncl_script_path}"' \
      'selMode="$selmode"' minVar=$minvar maxVar=$maxvar $MAP timeStep=$plotstep \
      'maskName="$mname"' $PROJSAT 'bStrg="$bstrg"' 'tStrg="$tstrg"' 'cdo="$CDO"'
eo_here

# run the ncl script
echo
#echo "=== Plotting ocean elevation ..."
echo "=== Plotting vertical velocity ..."
source ./scr_${EXP}_nclcmd.here
mkdir -p plots
echo "new Dir"
mv ${ofile}.eps plots/.
#check_error $? "In script icon_plot.ncl:"
#rm scr_${EXP}_nclcmd.here
echo "=== Done."

echo 
echo "=== Postprocessing finished for the hydrostatic ocean $EXP Shallow Water test."
echo "=== The plots can be found in "${dir}
