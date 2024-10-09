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

function diffWithExitCode {
#   set -x
  fileA=$1
  fileB=$2
  ofile=`mktemp`

  cdo -s diffv ${fileA} ${fileB} > $ofile

  nDiff=$(wc -l < $ofile)

  if (( $nDiff > 0 )); then
    echo "--------------------------------------------"
    echo "Files ${fileA} ${fileB} Differ!"
    cat $ofile;
    echo "--------------------------------------------"
  fi

  rm -f $ofile
  return $nDiff
}
function accumulationDiff {
  ifile=$1
  varName=$2
  maskName=$3

  diffWithExitCode "-div -selname,$varName $ifile -selname,$maskName $ifile " "-div -selname,${varName}_acc $ifile -selname,$maskName $ifile"
}

function directoryDiff {
#   set -x
  refDir=$1
  expDir=$2
  nDiff=0

  # abort if directory does not exist
  if [[ ! -d $refDir ]]; then
    return 99
  fi

  # abort if directory does not contain files to compare with to avoid false
  # positive results (0 == refCount after the diff)
  refCount=0

  # make sure that:
  # - only files (not directories) are used
  # - full path to the reference data is provided
  refList=$(find ${refDir} -maxdepth 1 -type f)

  for refFile in ${refList}; do 
    refFileBasename=$(basename ${refFile})
    case "${refFileBasename}" in
      *.nc*) # only netcdf files are taken into account
        DIFF='diffWithExitCode'
        ;;
      *)
        DIFF='true' # disabled for not nc-files
        ;;
    esac
    ${DIFF} ${refFile} ${expDir}/${refFileBasename}
    ((refCount = refCount + 1))
    if (( $nDiff > 0 )); then
      return $nDiff
    fi
  done

  if [[ 0 = ${refCount} ]]; then
    echo "Could not find anything to compare in ${refDir}!"
    return 99
  fi
}
