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

DATE=$1
builder=$2
builder_ID=$3

target_server=squall.zmaw.de
target_dir=/scratch/mpi/CC/mh0287/data/archive/${DATE}/buildbot/${builder}

ssh ${target_server} mkdir -p ${target_dir}
scp -r experiments ${target_server}:${target_dir}/${builder_ID}

