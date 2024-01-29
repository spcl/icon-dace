#! /bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#----------------------------------------------------------------------------------------------
usage()
{
    echo 
    echo usage: $0 options
    echo 
    echo This script fixes the ICON standard global attributes for DCMIP
    echo
    echo OPTIONS:
    echo   -h      show this message
    echo   -c      test case
    echo   -r      horizontal resolution
    echo   -l      levels
    echo   -g      grid
    echo   -t      time frequency
    echo   -d      description
    echo   -f      file to change
    echo 
}
#----------------------------------------------------------------------------------------------
while getopts ":c:r:l:g:t:d:f:h" opt
do
    case $opt in
	h) usage $0
	    exit 0
	    ;;
	c) case=$OPTARG
	    if [ "$case" == "" ]
	    then
		usage $0
	    fi
	    ;;
	r) hres=$OPTARG
	    if [ "$hres" == "" ]
	    then
		usage $0
	    fi
	    ;;
	l) levs=$OPTARG
	    if [ "levs" == "" ]
	    then
		usage $0
	    fi
	    ;;
	g) grid=$OPTARG
	    if [ "$grid" == "" ]
	    then
		usage $0
	    fi
	    ;;
	t) tfreq=$OPTARG
	    if [ "$tfreq" == "" ]
	    then
		usage $0
	    fi
	    ;;
	d) desc=$OPTARG
	    if [ "$desc" == "" ]
	    then
		usage $0
	    fi
	    ;;
	f) file=$OPTARG
	    if [ "$file" == "" ]
	    then
		usage $0
	    fi
	    ;;
	\?) echo "Invalid option: -$OPTARG" >&2
            exit 1
	    ;;
    esac
done
#----------------------------------------------------------------------------------------------
echo
if [ "$file" == "" ]
then
    echo ERROR: filename not given ...
    usage $0
    exit 1
fi
if [ "$case" == "" ]
then
    echo ERROR: case missing ...
    usage $0
    exit 1
fi
if [ "$hres" == "" ]
then
    echo ERROR: horizontal resolution missing ...
    usage $0
    exit 1
fi
if [ "$levs" == "" ]
then
    echo ERROR: levels missing ...
    usage $0
    exit 1
fi
if [ "$grid" == "" ]
then
    echo ERROR: grid missing ...
    usage $0
    exit 1
fi
if [ "$tfreq" == "" ]
then
    echo ERROR: time frequency missing ...
    usage $0
    exit 1
fi
if [ "$desc" == "" ]
then
    echo ERROR: description missing ...
    usage $0
    exit 1
fi
#----------------------------------------------------------------------------------------------
# first delete attributes, if they exist ...

ncatted -a Conventions,global,d,,"CF-1.0" $file 
ncatted -a model,global,d,, $file 
ncatted -a test_case,global,d,, $file 
ncatted -a horizontal_resolution,global,d,, $file 
ncatted -a levels,global,d,, $file 
ncatted -a grid,global,d,, $file 
ncatted -a time_frequency,global,d,, $file 
ncatted -a equation,global,d,, $file 
ncatted -a description,global,d,, $file 
#----------------------------------------------------------------------------------------------
# add the attributes
ncatted -a Conventions,global,c,c,"CF-1.0" $file 
ncatted -a model,global,c,c,"icon-mpi-dwd" $file
ncatted -a native_grid,global,c,c,"tri" $file  
ncatted -a test_case,global,c,c,"$case" $file 
ncatted -a horizontal_resolution,global,c,c,"$hres" $file 
ncatted -a levels,global,c,c,"$levs" $file 
ncatted -a grid,global,c,c,"$grid" $file 
ncatted -a time_frequency,global,c,c,"$tfreq" $file 
ncatted -a equation,global,c,c,"nonhydro" $file 
ncatted -a description,global,c,c,"$desc" $file 
#----------------------------------------------------------------------------------------------

