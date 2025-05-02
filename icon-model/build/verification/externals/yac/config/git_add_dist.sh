#!/bin/bash

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

set -eu
set -o pipefail

script_dir=$(cd "$(dirname "$0")"; pwd)
top_srcdir=$(cd "${script_dir}/.."; pwd)

cd "${top_srcdir}"

git update-index -q --refresh || {
  echo "ERROR: failed to update git index in '${top_srcdir}'" >&2
  exit 1
}

git diff-files --quiet || {
  echo "ERROR: '${top_srcdir}' has unstaged changes" >&2
  exit 1
}

git diff-index --cached --quiet HEAD -- || {
  echo "ERROR: '${top_srcdir}' has uncommited changes" >&2
  exit 1
}

./autogen.sh
./configure \
  --disable-mpi-checks \
  --enable-maintainer-mode \
  ac_cv_header_libfyaml_h=yes \
  ac_cv_header_mpi_h=yes \
  ac_cv_header_yaxt_h=yes \
  acx_cv_c_lib_func_MPI_Init= \
  acx_cv_c_lib_func_fy_node_mapping_lookup_scalar0_by_simple_key= \
  acx_cv_c_lib_func_xt_initialized= \
  acx_cv_fc_c_compatible_mpi=yes \
  acx_cv_fc_lib_func_MPI_INIT= \
  acx_cv_fc_module_MPI=yes \
  tj_cv_c_type_MPI_Fint=int \
  yac_cv_fc_is_contiguous_works=yes

distdir='yac-dist'
make distdir distdir="${distdir}"
for f in $(find "${distdir}" -type f -o -type l); do
  git add -f "${f#"${distdir}/"}";
done
rm -rf "${distdir}"
