#! /bin/sh

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

ser_vars2pass="my_iconppdir=\"${my_iconppdir}\";FC=\"${FC}\""
ser_funcsrc="$(readlink -f scripts/serialization/functions.sh)"
ser_setup="${ser_vars2pass};source ${ser_funcsrc}"

serialization_init() {
  echo "";
  echo "Preprocessing files for serialization ...";
  cd src/ || exit 1
  find ./ -type d -exec mkdir -p "../${my_iconppdir}/{}" \;
  echo "Updating SHA1 ..."
  find ./ -type f -regextype posix-egrep -regex ".*\.(f90|inc|incf|F90|cu|h)\$" -exec sh -c "${ser_setup};ser_update_file_hash \$1" _ {} \;
  echo "Updating headers ..."
  find ./ -type f \( -name "*.inc" -o -name "*.incf" \) -exec sh -c "${ser_setup};ser_update_inc_hash \$1" _ {} \;
  cd .. || exit 1
}

serialization_pp() {
  case $fortran_compiler in #(
    gcc*) :
      cd src/ || exit 1
      echo "Preprocessing files (CPP) ..."
      PP_FFLAGS="$(echo "${FFLAGS}" | sed -e 's|-I../pp ||' -e 's|-I../support ||' -e 's|-I../module ||')"
      gcc_setup="${ser_setup};PP_FFLAGS=\"${PP_FFLAGS}\""
      find ./ -type f \( -name "*.f90" -o -name "*.F90" \) -exec sh -c "${gcc_setup};ser_gcc_pp \$1" _ {} \;
      cd .. || exit 1
      ;; #(
    cray*) :
      local abs_include
      abs_include="$(readlink -f src/include)"
      cd src/ || exit 1
      echo "Preprocessing files (CPP) ..."
      PP_FFLAGS="$(echo "$FFLAGS" | sed -e 's/ -e Z //' -e 's/,ipa1//' -e 's|-I../pp ||' -e 's|-I../support ||' -e 's|-I../module ||')"
      cray_setup="${ser_setup};PP_FFLAGS=\"${PP_FFLAGS}\";abs_include=\"${abs_include}\""
      find ./ -type f \( -name "*.f90" -o -name "*.F90" \) -exec sh -c "${cray_setup};ser_cray_pp \$1" _ {} \;
      cd .. || exit 1
      ;; #(
    pgi*) :
      cd src/ || exit 1
      echo "Preprocessing files (CPP) ..."
      PP_FFLAGS="$(echo "$FFLAGS" | sed -e 's|-I../pp ||' -e 's|-I../support ||' -e 's|-I../module ||')"
      pgi_setup="${ser_setup};PP_FFLAGS=\"${PP_FFLAGS}\""
      find ./ -type f \( -name "*.f90" -o -name "*.F90" -o -name "*.cu" -o -name "*.h" \) -exec sh -c "${pgi_setup};ser_pgi_pp \$1" _ {} \;
      cd .. || exit 1
      ;; #(
    *) :
      ;;
  esac

  echo "Preprocessing files (Serialbox2) ...";
  sb2_setup="${ser_setup};SERIALBOX2_PPSER=\"${SERIALBOX2_PPSER}\""
  find "${my_iconppdir}" -type f \( -name "*.f90" -o -name "*.F90" \) -exec sh -c "${sb2_setup};ser_sb2_pp \$1" _ {} \;
}
