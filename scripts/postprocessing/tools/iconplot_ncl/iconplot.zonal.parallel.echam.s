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

#
#==========================================================================
# Driver script for postprocessing and visulization for ICON netcdf output
#==========================================================================
#
# History: 
# Taken from file HS_postpro_driver.bash:
#   Initial version by Hui Wan (MPI-M, 2010-07)
#   Pilar Ripodas (DWD, 2010-10), options for the Non Hydrostatic model
# iconplot.zonal.ncl-cdo.s:
#   Initial version by Martin Koehler (DWD, 2010-11)
#
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
#
# Parallel processing:
# (1) parallel computing:
#     pshell -p13 -f iconplot.zonal.parallel.list
#       which calls all variables like:
#       iconplot.zonal.parallel.s -v T
# (2) plotting:
#       iconplot.zonal.parallel.s
#==============================================================================
#
check_error()
{

# Check if the first parameter (return status) is not OK 

  if [ $1 -ne 0 ]; then

# Stop running this script and return the error status
    echo "ERROR: $2"
    exit $1
  fi
}

# argument list processing

while getopts v: option
do
   case $option in
    v) varname=$OPTARG;;    # variable name
   \?) errflg=1;;
   esac
done

#if [ $errflg -ne 0 ] ; then
#  echo "incorrect usage!" >&2
#  exit 1
#fi

#set -x


#==========================================================================
# quick setup (details below)

export expnum="exp12"
export expnum2="exp12"   #for comparison in plots, data must pre-exist
export Model="ICOHDC"    #ICONAM or ICOHDC
export EXPNAME="APE"     #APE or JWmoistnh

export top_title=${EXPNAME}"  "${expnum}"  ECHAM physics"
#export top_title=${EXPNAME}"  "${expnum}"  cosmo clouds, RRTM, ECHAMturb"
#export top_title=${EXPNAME}"  "${expnum}"  cosmo clouds, RRTM, ECHAMturb, SSO+orog"
#export top_title=${EXPNAME}"  "${expnum}"  cosmo clouds, Ritter, ECHAMturb"
#export top_title=${EXPNAME}"  "${expnum}"  no clouds, RRTM, Raschendorfer"
#export top_title=${EXPNAME}"  "${expnum}"   cosmo clouds, RRTM, Raschendorfer, Kmin=0.01m2/s"
#export top_title=${EXPNAME}"  "${expnum}"  clouds=1=grid-scale"
#export top_title=${EXPNAME}"  "${expnum}"  clouds=2=turbulence"
#export top_title=${EXPNAME}"  "${expnum}"  clouds=0=no cloud; no rad"
#export top_title=${EXPNAME}"  "${expnum}"  clouds=3=cosmo cloud; fixed"
#export top_title=${EXPNAME}"  "${expnum}"  clouds=4=new diagnostic clouds"
#export top_title=${EXPNAME}"  "${expnum}"  no cloud-radiation"

model_data_path="/e/uwork/mkoehler/icon/experiments/"${expnum}"/"
model_data_path2="/e/uwork/mkoehler/icon/experiments/"${expnum2}"/"
horizontal_resolution="R2B04"
vertical_resolution="L31"     # vertical levels
dtout=12                      # output frequency in [h]

time_stat="timselavg"         # time statistic: single time step or time period
noffset=400      # 180/540    # Number of input time steps skipped before the first time step range (CDO timavg)
nsets=400                     # Number of input time steps for output time step range               (CDO timavg)
time_stat_text="200-400 days" # DJF 90-180 days / JJA 270-360 days / 200-400 days
evol_istart=1                 # first file
evol_iend=40                  # last file

evol_stat="zonavg"            # evolution plot
#evol_stat="zonstd"

diag_climate=1
diag_evolution=1
diag_evolution_ps=0

if [ ! -z $varname ]; then
  do_computation=1            # parallel compuation(1=ON,0=OFF)
  do_computation_serial=0     # serial computation
  compute_remap_weights=0
  make_plot=0
else
  do_computation=0
  do_computation_serial=0     # optional: 1=ON,0=OFF (turn off when plotting only!)
  compute_remap_weights=0     # optional: 1=ON,0=OFF ( - " - )
  make_plot=1
fi
#==========================================================================



#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
# 0. About the dynamical core ( Hydrostatic or Non-Hydrostatic )
#--------------------------------------------------------------------------
# The variable Model can be ICOHDC (Hydrostatic) or ICONAM (Non-Hydrostatic)
# If set to ICOHDC the evolution of zonal temperature variance is 
# calculated and ploted at 750 hPa and the vertical variable in the plots 
# of the zonal mean climate is eta
# If set to ICONAM, the evolution of zonal temperature variance is 
# calculated and ploted at model level 24 and the vertical variable in the plots 
# of the zonal mean climate is the height
# 

#Model="ICOHDC"
#Model="ICONAM"

#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 What is the experiment ID? 
# The experiment ID is the first part of the file name of your model output, 
# which is followed by information on hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6) if the output has been split on several files.
# You may have used something like "A0001" just for simplicity.

#EXP="hat_heldsuarez"  
#EXP="JWmoistnh_DOM01"

# 1.2 In which directory are the data files located? (Don't forget the trailing "/"!)

#model_data_path="somewhere/"
#model_data_path="/e/uwork/mkoehler/icon/experiments/exp01/"
#model_data_path="/e/uscratch/mripodas/icon-dev/experiments/${EXPNAME}/"

# The ICON model developers have an automatic testing system which uses
# the following convention:
 
#dir=`pwd -P`
#icon_path=${dir%%scripts/ncl_scripts}
#model_data_path="${icon_path}experiments/$EXP/"

# 1.3 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31" 

config_string=${horizontal_resolution}${vertical_resolution}

# 1.4 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

export cell_type=3

# 1.5 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXPNAME} as part of the name
#  of the model output.)

#horizontal_resolution="R2B04"
#vertical_resolution="L31"        #ICOHDC
#vertical_resolution="L35"        #ICONAM  

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# A complete set of the contour plots includes
#
# a. evolution of zonal temperature variance at 750 hPa;
# b. the zonal mean "climate" computed for the user-specified period, including
#    > time and zonal mean temperature, u-wind, and v-wind,
#    > meridional flux of heat and zonal momentum caused by eddy activities,
#    > eddy kinetic energy.
# c. temporal evolution of the surface pressure glomal mean
#
# If you want to have data or plot of (a), set the following variable to 1

#diag_evolution=0

# If you want to have data or plot of (b), set the following variable to 1

#diag_climate=1

# If you want to have data or plot of (c), set the following variable to 1

#diag_evolution_ps=0

# Since the Held-Suarez test is an idealized climate simulation,
# this script assumes that the output has been split into a series of files. 
# The first file is named ${EXPNAME}"_RxBxxLxx_0001",the second one "_0002",
# and so on. Now specify the starting and ending file indices for plot (a)
# (see above)

#dtout=12        # output frequency in [h]
#evol_istart=1
#evol_iend=5

# Specify the starting and ending file indices for plot (b):

clim_istart=$evol_istart
clim_iend=$evol_iend

# Specify the starting and ending file indices for plot (c):

evol_ps_istart=$evol_istart
evol_ps_iend=$evol_iend

# Select time step to be plotted in zonal means (1 to nsteps)

#time_stat="timselavg" # time statistic: time average
#time_stat="timselvar" # time statistic: time variance
#time_stat_text="Time mean"
#time_stat_text="240 hours"

# If you want to obtain plots, set the following variable to 1 

#make_plot=1

# What format do you want, pdf, eps or ps?

export plot_file_format="ps"

# Before plotting, interpolation and statistics computation is needed. 
# If you have already done that and stored the data, turn this step off.

#do_computation=0   # (1=ON,0=OFF)

# To diagnose the zonal mean "climate", the ICON model output 
# on the geodesic grid needs to be interpolated to a Gaussian grid. 
# Choose the Gaussian grid by specifying the triangular truncation:

export trunc=85

# For the horizontal interpolation from the ICON grid to Gaussian grid,
# a large part of the time will be spent on calculating the remapping 
# weights. This may take very long when the resolution is high.
# Therefore we suggest calculating the weights only once and store them 
# for later use. 

#compute_remap_weights=0   # (1=ON,0=OFF)

# The remapping weights file generated by this script will be named, 
# e.g., icon_R2B04_spr0.90_tri_cell_to_T159.nc 
# If the weights are already available, specify the location:
# (Don't forget the trailing "/")

remap_weights_path="${model_data_path}remap_weights/"
#remap_weights_path='/e/uwork/mripodas/weights/'

# For the remapping, we also need to know which optimization was used 
# to generate the ICON grid. For example,
#   "ori"     : original icosahedral grid without optimization
#   "hro"     : Heikes-Randall optimization
#   "spr0.90" : spring dynamics, with spring coefficient 0.90

grid_optimization="spr0.90"

# Remove these files after finishing the diagnoses? 
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=0

# Where should the plot files be located? Don't forget the trailing "/".

plot_file_path="${model_data_path}plots/"

#--------------------------------------------------------------------------
# Now specify the directory in which the intermediate files 
# (excluding the remapping weights) should be placed. 
# Don't forget the trailing "/".

tmp_data_path="${model_data_path}tmp/"
tmp_data_path2="${model_data_path2}tmp/"

#--------------------------------------------------------------------------
# Do you want CDO to run in silence mode, or to report everything
# it is doing?

cdo_silence=1   #( 1 = silence mode; 0 = detailed report )

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

echo
echo "**********************************************************"
echo "***       ICON Tool Kit for Idealized Test Cases       ***"
echo "**********************************************************"
echo 

# Temporary variables

export script_path='/e/uhome/mkoehler/icon-dev/scripts/postprocessing/tools/'

#directory with the colomap files
export NCARG_COLORMAPS=${script_path}color_map:${NCARG_COLORMAPS}

resolution=${horizontal_resolution}${vertical_resolution}
fori=${model_data_path}${EXPNAME}_${resolution}
ftmp=${tmp_data_path}${EXPNAME}_${resolution}

if [ $cdo_silence -eq 1 ]; then
   silence='-s'
else
   silence=''
fi

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi


#==========================================================================
# Prepare remappping weights
#==========================================================================

weights=${remap_weights_path}"icon_"${horizontal_resolution}

if [ ${cell_type} -eq 3 ]; then
   weights=${weights}"_"${grid_optimization}"_cell_to_T"${trunc}".nc"
elif [ ${cell_type} -eq 6 ]; then
   weights=${weights}"_"${grid_optimization}"_vert_to_T"${trunc}".nc"
else
   echo "Wrong choice of cell_type. Should be 3 or 6 !"
   exit 1
fi

# compute weights

if [ ${compute_remap_weights} -eq 1 ]; then

 if [ ! -d ${remap_weights_path} ]; then
    mkdir -p ${remap_weights_path} 
 fi

 echo
 echo "=== Computing remapping weights (ICON to Gaussian) ..."

 label=$(printf "%04d" ${evol_istart})
 cdo $silence gendis,t${trunc}grid -seltimestep,1 -selname,PS ${fori}"_"${label}".nc" ${weights}

 check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'Computing remapping weights (ICON to Gaussian)'"
 echo "=== Done."

fi


#========================================================================
# Zonal mean climate: calculate and plot
#========================================================================

if [ $diag_climate -eq 1 ]; then

   if [ $do_computation -eq 1 ]; then

      echo
      echo "=== Computing statistics ..."

      #------------------
      # Select variables 
      #------------------

      ifile=$clim_istart
      while [ $ifile -le $clim_iend ] ; do

        label=$(printf "%04d" $ifile)
        #parallel# for varname in T U V W CC QV QC QI Q1 Q2 Q3 Q4 Q5 ; do
          cdo $silence selname,${varname}  ${fori}"_"$label".nc" ${ftmp}"_"$label"."${varname}".nc"
        #parallel# done

        check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'selname', file $label"
        ifile=` expr $ifile + 1 `
      done

      #-----------------------------------------
      # Merge all time steps into a single file
      #-----------------------------------------

      #parallel# for varname in T U V W CC QV QC QI Q1 Q2 Q3 Q4 Q5 ; do
          cdo $silence -r copy  ${ftmp}"_"[0-9][0-9][0-9][0-9]"."$varname".nc" \
                                ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc"
          rm ${ftmp}"_"[0-9][0-9][0-9][0-9]"."$varname".nc"
          check_error $? "In script iconplot.zonal.ncl-cdo.s: file merging"
      #parallel# done

      #-----------------------------------------
      # Compute time mean and variances (of all time steps) and select of a single time step
      #-----------------------------------------

      #parallel# for varname in T U V W CC QV QC QI Q1 Q2 Q3 Q4 Q5 ; do
        for opr in timselavg timselvar ; do
          cdo $silence $opr,$nsets,$noffset ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc"     \
                                            ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"."$opr".nc"
        done
      #parallel# done

      check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'mean and variance'"

   fi

   if [ $do_computation_serial -eq 1 ]; then

      #-----------------------------------------
      # Eddy kinetic energy
      #-----------------------------------------

      cdo $silence mulc,0.5 -add ${ftmp}"_"$clim_istart"-"$clim_iend".U.timselvar.nc" \
                                 ${ftmp}"_"$clim_istart"-"$clim_iend".V.timselvar.nc" \
                                 ${ftmp}"_"$clim_istart"-"$clim_iend".eddy-KE.nc"

      rm ${ftmp}"_"$clim_istart"-"$clim_iend".U.timselvar.nc"
      rm ${ftmp}"_"$clim_istart"-"$clim_iend".V.timselvar.nc"

      check_error $? "In script iconplot.zonal.ncl-cdo.s: eddy kinetic energy"

      #-----------------------------------------
      # Eddy momentum/heat flux
      #-----------------------------------------

      for varname in T U ; do

          cdo $silence timselavg,$nsets,$noffset -mul ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc" \
                                                      ${ftmp}"_"$clim_istart"-"$clim_iend".V.nc"          \
                                                      ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timselavg.nc"

          cdo $silence mul ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".timselavg.nc" \
                           ${ftmp}"_"$clim_istart"-"$clim_iend".V.timselavg.nc"          \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"

          cdo $silence sub ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timselavg.nc" \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"      \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"pVp.nc"

          rm ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timselavg.nc"
          rm ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"

          check_error $? "In script iconplot.zonal.ncl-cdo.s: eddy flux of $varname"
      done

      #-------------------------------------------------
      # Interpolate to Gaus grid and compute zonal mean
      #-------------------------------------------------

      varlabel1=("U."$time_stat  "V."$time_stat  "OMEGA."$time_stat "T."$time_stat    \
                 "Qv."$time_stat "Qw."$time_stat "Qi."$time_stat    "ACLC."$time_stat \
                 "ACLCAC."$time_stat             "Q4."$time_stat    "Q5."$time_stat   \
                 "UpVp"          "TpVp"          "T.timselvar"      "eddy-KE"         )
      varlabel2=("U"             "V"             "OMEGA"            "T"               \
                 "Qv"            "Qw"            "Qi"               "ACLC"            \
                 "ACLCAC"        "Q4"            "Q5"                                 \
                 "UpVp"          "TpVp"          "TpTp"             "eddy-KE"         )

      nvar=${#varlabel1[*]} 
      ivar=0

      while [ $ivar -lt $nvar ]; do

        cdo $silence zonavg -remap,t${trunc}grid,${weights}    \
          ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel1[$ivar]}".nc" \
          ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel2[$ivar]}".zmta.T"$trunc".nc"

        check_error $? "In script iconplot.zonal.ncl-cdo.s: remap and zonal mean"
        rm ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel1[$ivar]}".nc"


      ivar=` expr $ivar + 1 `
      done
   fi

   if [ $make_plot -eq 1 ]; then

      export DataPath=${tmp_data_path}
      export DataPath2=${tmp_data_path2}
      export DataID=${EXPNAME}_${resolution}"_"${clim_istart}-${clim_iend}
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export Resolution=${horizontal_resolution}${vertical_resolution}
      export trunc
      export dtout
      export time_stat_text

      ncl iconplot.zonal-vert.echam.ncl
   fi
fi


#==================================================================================================
# Time series of zonal mean or variance at single levels: calculate and plot.
#==================================================================================================

if [ ${diag_evolution} -eq 1 ]; then

  if [ ${do_computation} -eq 1 ]; then

   #for varname in QV ; do
    #parallel# for varname in T U V W CC QV QC QI Q1 Q2 Q3 Q4 Q5 ; do
   #for level in 80 69 52 38 11 ; do    # synchonize with iconplot.zonal-time.ncl
    for level in 30 20 10 5 ; do    # synchonize with iconplot.zonal-time.ncl

      if [ $Model == "ICOHDC" ]; then
       #templ="temp750"
        templ=$varname"_L"$level
      fi
      if [ $Model == "ICONAM" ]; then
        templ=$varname"_L"$level
      fi

      #-------------------------
      # level and variable selection (vertical interpolation)
      #-------------------------
   
      ifile=$evol_istart
      while [ $ifile -le $evol_iend ] ; do
   
        label=$(printf "%04d" $ifile)
        if [ $Model == "ICOHDC" ]; then
        #cdo $silence ml2pl,75000 -selname,T,PS,PHIS ${fori}"_"$label".nc" ${ftmp}"_"$label"."${templ}".nc"
         cdo $silence sellevel,$level -selname,$varname ${fori}"_"$label".nc" ${ftmp}"_"$label"."${templ}".nc"
         check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'Interpolating model temperature to 750 hPa', file $label"
        fi
        if [ $Model == "ICONAM" ]; then
         cdo $silence sellevel,$level -selname,$varname ${fori}"_"$label".nc" ${ftmp}"_"$label"."${templ}".nc"
         check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'Selecting variable at a model level', file $label"
        fi
        ifile=` expr $ifile + 1 `
      done
   
      #--------------------------
      # merge to a single file 
      #--------------------------
   
      cdo $silence -r copy ${ftmp}"_"[0-9][0-9][0-9][0-9]"."${templ}".nc" \
                           ${ftmp}"_"${evol_istart}-${evol_iend}"."${templ}".nc"
   
      check_error $? "In script iconplot.zonal.ncl-cdo.s: file merging"
   
      rm ${ftmp}"_"[0-9][0-9][0-9][0-9]"."${templ}".nc"
   
      #-----------------------------
      # interpolation to Gaus grid
      #-----------------------------
   
      echo
      echo "=== Computing "$evol_stat" of "$varname" at level "${level}"..."
   
      cdo $silence -r -$evol_stat -remap,t${trunc}grid,${weights} \
          ${ftmp}"_"${evol_istart}-${evol_iend}"."${templ}".nc" \
          ${ftmp}"_"${evol_istart}-${evol_iend}".T"${trunc}"."${templ}".zonvar.nc"
   
      check_error $? "In script iconplot.zonal.ncl-cdo.s: part 'Computing zonal variance of temperature at '"${level}

    done
    #parallel# done
  fi
  
  if [ ${make_plot} -eq 1 ]; then

      export DataPath=${tmp_data_path}
      export DataPath2=${tmp_data_path2}
      export DataID=${EXPNAME}_${resolution}"_"${evol_istart}-${evol_iend}".T"${trunc}
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export Resolution=${horizontal_resolution}${vertical_resolution}
      export dtout
      export evol_stat
      export trunc
      export level

      echo "=== Plotting temperature variance ..."
      ncl iconplot.zonal-time.echam.ncl
      check_error $? "In scripte iconplot.zonal.ncl-cdo.s: part 'Plotting temperature variance'"
  fi

fi


#========================================================================
# Surface pressure global mean: calculate and plot.
#========================================================================

if [ ${diag_evolution_ps} -eq 1 ]; then

   if [ ${do_computation_serial} -eq 1 ]; then

      #----------------------------------
      # select PS and calculate global mean
      #----------------------------------

      ifile=$evol_ps_istart
      while [ $ifile -le $evol_ps_iend ] ; do

        label=$(printf "%04d" $ifile)

        cdo $silence fldmean  -selname,PS ${fori}"_"$label".nc" ${ftmp}"_"$label".PSmean.nc"       

        ifile=` expr $ifile + 1 `

      done
   
      #--------------------------
      # merge to a single file 
      #--------------------------

      cdo $silence -r copy ${ftmp}"_"[0-9][0-9][0-9][0-9]".PSmean.nc" \
                           ${ftmp}"_"${evol_ps_istart}-${evol_ps_iend}".PSmean.nc"
   
      check_error $? "In script iconplot.zonal.ncl-cdo.s: file merging"
   
      rm ${ftmp}"_"[0-9][0-9][0-9][0-9]".PSmean.nc"
   fi

   if [ ${make_plot} -eq 1 ]; then

      export DataPath=${tmp_data_path}
      export DataID=${EXPNAME}_${resolution}"_"${evol_ps_istart}-${evol_ps_iend}".PSmean"
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export VarName="PS"
      export Resolution=${horizontal_resolution}${vertical_resolution}

      echo "=== Plotting Surface Pressure evolution ..."
      ncl ${script_path}HS_plot_PS-mean_time_series.ncl 
      check_error $? "In scripte iconplot.zonal.ncl-cdo.s: part 'Plotting Surface Pressure global mean'"
      echo "=== Done."
   fi

fi


#========================================================================
# Clean up
#========================================================================

echo "=== The plots can be found in "${plot_file_path}

if [ $rm_tmp_files -eq 1 ]; then

   rm ${tmp_data_path}/${EXPNAME}_${resolution}*.nc 

   if [ `ls ${tmp_data_path} |wc -l` -eq 0 ]; then
      rm -rf ${tmp_data_path}
   fi 
   echo "=== The temporary data have been removed."

else
   echo "=== The pressure level data and spectral coefficients can be found in "${tmp_data_path}
fi

exit
