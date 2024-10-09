#!/bin/bash

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

#-----------------------------------------------------------------------------
# Script to check out the scripts for ICON-Land initial files generation from
# the ICON-Land repository at gitlab.dkrz.de:jsbach/jsbach . Scripts are put
# into the subdirectory "scripts" which must not exist.

set -e

branch="land4icon-mpim" # branch from which to use the scripts

if [[ -d scripts ]]; then
  echo "$0: subdirectory 'scripts' already exists ... exiting"
  exit 1
fi

git clone -n -b master git@gitlab.dkrz.de:jsbach/jsbach scripts >& /dev/null
#git clone -n -b master https://gitlab.dkrz.de/jsbach/jsbach.git scripts >& /dev/null
cd scripts
git config core.sparseCheckout true >& /dev/null
echo /scripts/preprocessing/initial_files >> .git/info/sparse-checkout
git checkout ${branch} >& /dev/null
mv scripts/preprocessing/initial_files/* ./
rm -rf scripts
cd ..

echo ""
echo ">>> Checked out the scripts for ICON-Land initial file generation from branch ${branch}".
echo "    cd to ./scripts and see README_initial_files.md for further instructions."
echo ""

