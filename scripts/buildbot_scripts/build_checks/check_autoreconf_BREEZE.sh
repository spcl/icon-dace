#/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# Checks that the configure scipt is generated with the vanilla Autoconf 2.69
# and matches configure.ac.
# Must be run on Breeze where we have the vanilla Autoconf 2.69.

set -eu
set -o pipefail

icon_dir=$(unset CDPATH; cd "$(dirname "$0")/../../.."; pwd)
tmp_dir=$(mktemp -d)
trap 'rm -rf -- "${tmp_dir}"' EXIT

# Make a copy of the configure script in the repo:
cp ${icon_dir}/configure ${tmp_dir}/configure

# Re-generate the configure script:
(
  case $(lsb_release -sc) in
    stretch)
      sw_root='/data/mpi/sclab/sip/m300488/sw/gcc-6.3.0'
      autoconf_root="${sw_root}/autoconf-2.69-6s27zjm"
      # We need Automake for the aclocal script:
      automake_root="${sw_root}/automake-1.16.1-vcfcb2c"
      export PATH="${autoconf_root}/bin:${automake_root}/bin:$PATH"
      ;;
    bullseye)
      module purge
      module load autoconf/2.69 automake
      ;;
  esac
  cd ${icon_dir}
  autoreconf -fvi
)

diff -u "${tmp_dir}/configure" "${icon_dir}/configure" >&2 || {
  cat >&2 <<_EOF
ERROR: the configure script is inconsistent with the contents of configure.ac
#      (see above)
_EOF
  exit 1
}
