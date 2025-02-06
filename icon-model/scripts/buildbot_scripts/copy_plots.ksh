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

# This script is used to copy the eps plot files to mpipc91.
# On mpipc91 the eps-files are converted to jpg-file.
#-------------------------------------------------------------------
#
# First version by Walter Sauf (MPI-M, 2010-02-26)
#
# $Rev$:     Revision of last commit
# $Author$:  Author of last commit
# $Date$:    Date of last commit
#-------------------------------------------------------------------
set -x
REV=$1
SLAVE=$2
BUILDER=$3
BUILDER_NR=$4
#DATE=`date '+%Y-%m-%d'`
W_DIR=$5
DATE=$6
BRANCH=$7
PUPDATE=`date '+%d.%m %H:%M'`


#===============================================================
# functions
#===============================================================

copy_files ()
{

  for file in $FILES
  do
    dir_name=`echo $file | cut -d '/' -f2-10`
    name=`basename $dir_name .$1`
    dir=`dirname $dir_name`
    mkdir -p ${BASE_DIR}/${dir}

 
    if [ "x$1" = "xps" -o "x$1" = "xeps" ]
    then
      if [ ${CONVERT_STATUS} = 0 ]
      then
        ${CONVERT} experiments/${dir}/${name}.$1 ${BASE_DIR}/${dir}/${name}.png
        echo "INFO: Convert ${name}.$1  --> ${name}.png "
      else
        echo "INFO: No convert available " 
      fi
    else
      cp experiments/${dir}/${name}.$1 ${BASE_DIR}/${dir}/${name}.$1
    fi
  done
}

build_png_files ()
{

  for file in $FILES
  do
    dir_name=`echo $file | cut -d '/' -f2-10`
    name=`basename $dir_name .$1`
    dir=`dirname $dir_name`

    mkdir -p ${BASE_DIR}/${dir}

    if [ ${CONVERT_STATUS} = 0 ]
    then
      ${CONVERT} experiments/${dir}/${name}.$1 ${BASE_DIR}/${dir}/${name}.png
      echo "INFO: Convert ${name}.$1  --> ${name}.png "
    else
      cp experiments/${dir}/${name}.$1 ${BASE_DIR}/${dir}/${name}.$1
      echo "INFO: Copy ${name}.$1"
    fi
  done
}

#===============================================================
# Script begin
#===============================================================
#set -x
#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Save Plots for Buildbot web-Page
#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Check if convert is known

CONVERT=`which convert`
CONVERT_STATUS=$?

if [ "x${BUILDER}" = "x" ]
then
  echo "!!!! BUILDER not set "
  exit 1
fi
#===================================================================
#
# This part creates the directory for the old buildbot-plot webpages
#
#===================================================================

#if [ -d /tmp/${BUILDER} ]
#then
#  rm -rf /tmp/${BUILDER} 
#fi

#echo "BB_SYSTEM=${SLAVE}"
#echo "BB_SLAVE=${BUILDER}"

#REV=`svn info | grep Revision | cut -d ':' -f2`
#echo "REV=${REV}"

#mkdir /tmp/${BUILDER} 
#echo "_COMPUTER_ ${SLAVE}" > /tmp/${BUILDER}/job_info.txt
#echo "_BUILDER_ ${BUILDER}" >> /tmp/${BUILDER}/job_info.txt
#echo "_REVISION_ ${REV}" >> /tmp/${BUILDER}/job_info.txt
#echo "_UPDATE_ ${PUPDATE}" >> /tmp/${BUILDER}/job_info.txt
#echo "_BUILDID_ ${BUILDER_NR}" >> /tmp/${BUILDER}/job_info.txt
#pwd
#find experiments -name '*.eps' -exec cp {} /tmp/${BUILDER}/. \;

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Save Plots for archive
#---------------------------------------------------------------------
#---------------------------------------------------------------------

BASE_DIR=/tmp/BuildBot/${BUILDER}/archive/${DATE}/buildbot/${REV}/${BRANCH}/${BUILDER}/${BUILDER_NR}

# This directory is only a temporary directory. Buildbot is downloadin the directpory
# to mpipc91 and saves the data in /pool/data/ICON/archiv/${DATE}/buildbot/${REV}/${BUILDER}/${BUILDER_NR}
# If a new run is started the BASE_DIR can be deleted. In the other case we would always create new directories in /mp


if [ -d  /tmp/BuildBot/${BUILDER}/archive ]
then
  rm -rf /tmp/BuildBot/${BUILDER}/archive
fi

mkdir -p ${BASE_DIR}

pwd
echo "Copy eps-files"
FILES=`find experiments -name '*.eps'`
copy_files eps

pwd
echo "Copy ps-files"
FILES=`find experiments -name '*.ps'`
copy_files ps

echo "Copy png-files"
FILES=`find experiments -name '*.png'`
copy_files png

#==================== End =====================================================

