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

script_dir=$(cd "$(dirname "$0")"; pwd)

in_file="${script_dir}/icon-config-doc-depgraph.dot"

out_format='svg'
out_file="${script_dir}/icon-config-doc-depgraph.${out_format}"

dot -Gconcentrate=true -T"${out_format}" "${in_file}" -o "${out_file}"
