# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#=======================================================================
# Variable registry for the ICON APE postprocessing package
#
# Author: Hui Wan (MPI, 2011-04)
# Author: Marc Salzmann (MPI, 2011-04)
#
# This scripts builds up a registry (table) containing information
# needed for data processing (include variable name, GRIB code, 
# unit, etc.) and for plotting (e.g., contour levels, color map). 
#
# When user calls this script for a field by giving its name in 
# the ICON model output, e.g., like
#  . ./lookup_variable.ksh PS
# the script searches through the registry and returns the entry index 
# of the field as variable "ie" to the calling script. If no matching
# entry is found, "ie" gets the value -1.
#=======================================================================

ie=-1   # initialize the result to be returned
je=-1   # loop index, internal use only

#---------------------------------------
# Build registry
#---------------------------------------

je=`expr $je + 1 `
shortname[je]="PS"
longname[je]="Surface pressure (Pa)"
varcode[je]=134
vartble[je]=128
tmpname[je]="aps"
afterbn[je]=0
plotscale[je]=0.01
plotmin[je]=940
plotmax[je]=1040
plotint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffmax[je]=50
diffint[je]=5
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="PHIS"
longname[je]="Surface geopotential (m2/s2)"
varcode[je]=129
vartble[je]=128
tmpname[je]="geosp"
afterbn[je]=0
plotscale[je]=1
plotmin[je]=-10
plotmax[je]=10
plotint[je]=1
diffmax[je]=10
diffint[je]=1
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="T"
longname[je]="Temperature (K)"
varcode[je]=130
vartble[je]=128
tmpname[je]="st"
afterbn[je]=1
plotscale[je]=1
plotmin[je]=160
plotmax[je]=305
plotint[je]=5
diffmax[je]=20
diffint[je]=2
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="U"
longname[je]="Zonal velocity (m/s)"
varcode[je]=131
vartble[je]=128
tmpname[je]="var131"
afterbn[je]=1
plotscale[je]=1
plotmin[je]=-20
plotmax[je]=120
plotint[je]=10
diffmax[je]=40
diffint[je]=4
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="V"
longname[je]="Meridional velocity (m/s)"
varcode[je]=132
vartble[je]=128
tmpname[je]="var132"
afterbn[je]=1
plotscale[je]=1
plotmin[je]=-5
plotmax[je]=5
plotint[je]=0.5
diffmax[je]=1.5
diffint[je]=0.3
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="OMEGA"
longname[je]="Vertical velocity (Pa/s)"
varcode[je]=135
vartble[je]=128
tmpname[je]="var135"
afterbn[je]=1
plotscale[je]=1
plotmin[je]=-0.3
plotmax[je]=0.3
plotint[je]=0.03
diffmax[je]=0.05
diffint[je]=0.01
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Qv"
longname[je]="Specific density of water vapour (g/kg)"
varcode[je]=133
vartble[je]=128
tmpname[je]="q"
afterbn[je]=0
plotscale[je]=1e3
plotmin[je]=1
plotmax[je]=15
plotint[je]=2
diffmax[je]=0.5
diffint[je]=0.1
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="QV"
longname[je]="Specific density of water vapour (g/kg), diagnostic"
varcode[je]=51
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e3
plotmin[je]=1
plotmax[je]=15
plotint[je]=2
diffmax[je]=0.5
diffint[je]=0.1
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Q1"
longname[je]="Specific density of water vapour (g/kg), grid-scale"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e3
plotmin[je]=1
plotmax[je]=15
plotint[je]=2
diffmax[je]=0.5
diffint[je]=0.1
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Qw"
longname[je]="Specific density of cloud water (kg/kg)"
varcode[je]=153
vartble[je]=128
tmpname[je]="xl"
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=10
plotmax[je]=190
plotint[je]=20
diffmax[je]=50
diffint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="QC"
longname[je]="Specific density of cloud water (mg/kg), diagnostic"
varcode[je]=31
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=10
plotmax[je]=190
plotint[je]=20
diffmax[je]=50
diffint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Q2"
longname[je]="Specific density of cloud water (mg/kg), grid-scale"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=10
plotmax[je]=190
plotint[je]=20
diffmax[je]=50
diffint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Qi"
longname[je]="Specific density of cloud ice (kg/kg)"
varcode[je]=154
vartble[je]=128
tmpname[je]="xi"
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=2
plotmax[je]=100
plotint[je]=5
diffmax[je]=20
diffint[je]=2
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="QI"
longname[je]="Specific density of cloud ice (mg/kg), diagnostic"
varcode[je]=33
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=10
plotmax[je]=190
plotint[je]=20
diffmax[je]=50
diffint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Q3"
longname[je]="Specific density of cloud ice (mg/kg), grid-scale"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=10
plotmax[je]=190
plotint[je]=20
diffmax[je]=50
diffint[je]=10
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Q4"
longname[je]="Specific density of rain (mg/kg), grid-scale"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=2
plotmax[je]=100
plotint[je]=5
diffmax[je]=20
diffint[je]=2
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="Q5"
longname[je]="Specific density of snow (mg/kg), grid-scale"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e6
plotmin[je]=5
plotmax[je]=200
plotint[je]=10
diffmax[je]=40
diffint[je]=4
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="ACLC"
longname[je]="Cloud cover (%)"
varcode[je]=162
vartble[je]=128
tmpname[je]="aclc"
afterbn[je]=0
plotscale[je]=1e2
plotmin[je]=5
plotmax[je]=80
plotint[je]=10
diffmax[je]=30
diffint[je]=3
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="CC"
longname[je]="Cloud cover [%]"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1e2
plotmin[je]=5
plotmax[je]=80
plotint[je]=10
diffmax[je]=30
diffint[je]=3
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]="W"
longname[je]="Vertical velocity (m/s)"
varcode[je]=999
vartble[je]=201
tmpname[je]=""
afterbn[je]=0
plotscale[je]=1
plotmin[je]=-0.03
plotmax[je]=0.03
plotint[je]=0.003
diffmax[je]=0.005
diffint[je]=0.001
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

je=`expr $je + 1 `
shortname[je]=""
unit[je]=""
varcode[je]=
vartble[je]=
tmpname[je]=""
afterbn[je]=
plotscale[je]=
plotmin[je]=
plotmax[je]=
plotint[je]=
diffmax[je]=
diffint[je]=
colormap[je]="WhBlGrYeRe"
colorstart[je]=0
colorend[je]=101
diffcolormap[je]="ViBlGrWhYeOrRe"
diffcolorstart[je]=2
diffcolorend[je]=250

#-----------------------
nentry=${#shortname[*]}
#echo Defined a look-up table of $nentry entries.



#-------------------------------------------------------------------
# Search through the registry for the user-specified variable name
#-------------------------------------------------------------------

je=0
until [[ $je -eq $nentry ]]; do
  if [ "$1" == "${shortname[$je]}" ]; then
    echo Found matching entry \#$je \(${shortname[$je]}\) in the registry.
    ie=$je
    break
  fi
  je=`expr $je + 1 `
done

if [ $ie -eq -1 ]; then
  echo Did not find variable $1 in the table!
  exit
fi

