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

#-------------------------------------------------------------------
# Driver for interactive visualization: the main idea is to exploit
# the shell functionalities to handle the user options, including
# default values, and then pass all the information to the ncl scripts
# by setting the required environmental variables.
#
# Use the flag -h to see the documentation.
#
# First version by Marco Restelli (MPI-M, 2009-03-31)
#-------------------------------------------------------------------

version="0.0.1"

# Here it is possible to set the defaults --------------------------

# terminal
default_output_terminal="x11"

# time level
default_time_level=-1

# fixed dimension
default_fix_dimension="lat"

# fixed level in above dimension
default_fix_level=1

# colormap
default_colormap="testcmap"
#"BlWhRe"

# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Initialize the defaults where possible
output_terminal=${default_output_terminal}
time_level=${default_time_level}
fix_level=${default_fix_level}
fix_dimension=${default_fix_dimension}
colormap=${default_colormap}
# ------------------------------------------------------------------


while getopts "C:FL:D:S:T:ho:t:vz:" Option
 do
  case $Option in
   h ) echo
            #------------------------------------------------------------------
       echo "Usage: `basename $0`"\
        "[-h] [-C colormap] [-F] [-L start:end:ncontours] -S scalar_varname [-T] [-o output_file] [-t time_lev] [-v] [-D dimension] [-z level] input_file.nc"
       echo
       echo "Plot data from the NetCDF file input_file.nc using the NCAR"
       echo "command language interface. Edit this file to change the defaults."
       echo
       echo "Command line flags (flags in square braces are optional):"
       echo "  -h  show this help and exit"
       echo "  -C  colormap; default: \"${default_colormap}\";"
       echo "      some other possibilities are: \"hotres\" \"testcmap\""
       echo "      \"gsdtol\" \"helix\""
       echo "  -F  force overwriting of existing files"
       echo "  -L  contour levels in the form start:end:ncontours"
       echo "  -S  specify the name of the scalar variable [MANDATORY!]"
       echo "  -T  terminal; possible values \"x11\" \"ncgm\" \"ps\" \"eps\"" 
       echo "      \"epsi\" \"pdf\""
       echo "      default: \"${default_output_terminal}\""
       echo "  -o  output_file; default constructed from the input file"
       echo "  -t  time level; this can be specified in three ways:"
       echo "      -t 5      ->  time step 5"
       echo "      -t 5:8    ->  time steps from 5 to 8, increment 1"
       echo "      -t 5:-1:2 ->  time steps from 5 to end, increment 2"
       echo "      use -1 as start or end index to indicate the last time"
       echo "      step in the input file; values start from 0;"
       echo "      when specifing a range, the output terminal must support"
       echo "      multiple pages; default: ${default_time_level}"
       echo "  -v  set verbose (use this also to see errors from ncl)"
       echo "  -D  fixed dimension of the section; possible value \"lat\" \"lon\" \"lev\""
       echo "      default: \"${default_fix_dimension}\""
       echo "  -z  axis level; default: ${default_fix_level}"
       echo "      (values starting from 1)"
       echo
       echo "REMARK: visualization of vector fields and superposition of"
       echo "multiple scalar fileds not yet implemented."
       echo
       exit;;
   C ) have_colormap=1
       colormap="${OPTARG}";;
   F ) FORCE=1;;
   L ) have_contour_levels=1
       contour_levels="${OPTARG}";;
   D ) have_fix_dimension=1
       fix_dimension="${OPTARG}";;
   S ) have_varname=1
       varname="${OPTARG}";;
   T ) have_output_terminal=1
       output_terminal="${OPTARG}";;
   o ) have_output_file=1
       output_file="${OPTARG}";;
   t ) have_time_level=1
       time_level="${OPTARG}";;
   v ) verbose=1;;
   z ) have_fix_level=1
       fix_level="${OPTARG}";;
  esac
done
shift $(($OPTIND - 1))

# ------------------------------------------------------------------
# Function definitions

function get_numeric_range()
{
# Given a string in the form  aa:bb:cc  extract the tree numbers aa,
# bb and cc. The string can include white spaces, and each number can
# have a sign, decimals and an exponent. Valid numbers are, for
# instance, "8", "+8", "-8", "-8.", "-8.54", "-8.54e2", "-8.54E-6"
#
# The output are the two variables range range_n.

 local wrange my_numeric_format sed_extract n sed_extracted

 wrange=$1 # working variable

 # 0) define some search patterns
 my_numeric_format='[+-]\{0,1\}[0-9]*\.\{0,1\}[0-9]*[eE]\{0,1\}[+-]\{0,1\}[0-9]*'
 sed_extract="s/\(^${my_numeric_format}\)\(.*\)/\1/"

 # 1) strip spaces
 wrange=`echo ${wrange} | sed "s/ //g"`

 # 2) extract the three numbers
 n=0
 sed_extracted="x"
 until [ "${sed_extracted}" = "" ]
  do
   sed_extracted=`echo ${wrange} | sed "${sed_extract}"`
   wrange=`echo ${wrange} | sed "s/^${sed_extracted}:\{0,1\}//"`
   if [ ! "${sed_extracted}" = "" ]
    then
     n=$((n+1))
     range[${n}]=${sed_extracted}
   fi
 done   
 range_n=${n}
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Check input and output files

# Check input file
case $# in
 "1" );; # everything is fine
  *  ) echo "Wrong number of input files: $#"
       exit 1;;
esac
if [ ! -e $1 ]
 then
  echo "Input file $1 does not exist!"
  exit 1
fi
input_file=$1

# Define the output file
if [ ! ${have_output_file} ]
 then
  file_suffix="_itime${time_level}_var${varname}_zlev${vert_level}"
  # the extension is added by ncl
  output_file="`echo ${input_file} | sed "s/.nc$/${file_suffix}/"`"
  # ncl does not like ":" or " " in the output file
  output_file="`echo ${output_file} | sed "s/:/-/g"`"
  output_file="`echo ${output_file} | sed "s/ //g"`"
fi
# Check the output file
if [ -e ${output_file}.${output_terminal} ]
 then
  if [ ! ${FORCE} ]
   then
    echo "Output file ${output_file}.${output_terminal} already exists!"
    echo "Use -F to overwrite existing files"
    exit 1
  fi
fi
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Check variable name
if [ ! ${have_varname} ]
 then
  echo "Please specify the variable name by -Svarname!"
  exit 1
fi
# add here a test to check whether the variable appears in the input file
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Check contour levels
contour_levels_ncont=-1
if [ ${have_contour_levels} ]
 then
  get_numeric_range "${contour_levels}"
  contour_levels_start=${range[1]}
  contour_levels_end=${range[2]}
  contour_levels_ncont=${range[3]}
fi
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Check time levels
get_numeric_range "${time_level}"
time_level_start=${range[1]}
case ${range_n} in
 "1" ) time_level_end=${time_level_start}
       time_level_step="1";;
 "2" ) time_level_end=${range[2]}
       time_level_step="1";;
 "3" ) time_level_end=${range[2]}
       time_level_step=${range[3]};;
  *  ) echo "Wrong time level specification: ${time_level}"
       exit 1;;
esac
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Start the ncl script
if [ ${verbose} ]
 then
  echo "Starting the ncl script \"plot_single_contour.ncl\" with variables:"
  echo "  varname         = ${varname}"
  echo "  input_file      = ${input_file}"
  echo "  time_level_start= ${time_level_start}"
  echo "  time_level_end  = ${time_level_end}"
  echo "  time_level_step = ${time_level_step}"
  echo "  fix_level       = ${fix_level}"
  echo "  fix_dimension   = ${fix_dimension}"
  echo "  colormap        = ${colormap}"
  echo "  contour_levels_start = ${contour_levels_start}"
  echo "  contour_levels_end   = ${contour_levels_end}"
  echo "  contour_levels_ncont = ${contour_levels_ncont}"
  echo "  output_terminal = ${output_terminal}"
  echo "  output_file     = ${output_file}"
fi

export varname
export input_file
export time_level_start
export time_level_end
export time_level_step
export fix_level
export fix_dimension
export colormap
export contour_levels_start
export contour_levels_end
export contour_levels_ncont
export output_terminal
export output_file

#ncl_output=$(ncl `dirname $0`/plot_single_contour.ncl)
ncl `dirname $0`/plot_torus_section.ncl

if [ ${verbose} ]
 then
  echo "${ncl_output}"
fi
# ------------------------------------------------------------------

