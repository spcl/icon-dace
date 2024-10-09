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

ser_update_file_hash() {
  local hash_file="../${my_iconppdir}/$1.sha1"
  local modif_file="../${my_iconppdir}/$1.modified"
  if [ -e "${hash_file}" ] && ( sha1sum --status -c "${hash_file}" );then
    rm -f "${modif_file}" 
  else
    rm -f "${hash_file}";
    sha1sum "$1" > "${modif_file}"
  fi
}

ser_update_inc_hash() {
  local pp_file="../${my_iconppdir}/$1"
  local hash_file="../${my_iconppdir}/$1.sha1"
  local modif_file="../${my_iconppdir}/$1.modified"
  if [ -e "${modif_file}" ];then
    cp "$1" "${pp_file}"
    mv "${modif_file}" "${hash_file}"
  fi
}

ser_gcc_pp() {
  local pp_file="../${my_iconppdir}/$1"
  local hash_file="../${my_iconppdir}/$1.sha1"
  local modif_file="../${my_iconppdir}/$1.modified"
  if [ -e "${modif_file}" ];then
    if grep -q '!$ser' $1; then
      ${FC} -I./include ${PP_FFLAGS} -o "${pp_file}" -C -E "$1"
    else
      cp "$1" "${pp_file}"
      mv "${modif_file}" "${hash_file}"
    fi
  fi
}

ser_cray_pp() {
  local pp_file="../${my_iconppdir}/$1"
  local hash_file="../${my_iconppdir}/$1.sha1"
  local modif_file="../${my_iconppdir}/$1.modified"
  if [ -e "${modif_file}" ];then
    if grep -q '!$ser' $1; then
      cur_pwd="$(pwd)"
      cd "$(dirname "$1")" || exit 1
      ${FC} -I"${abs_include}" ${PP_FFLAGS} -eP "$(basename "$1")"
      cd "${cur_pwd}" || exit 1
      mv "$(echo "$1" | sed -e 's/.f90//I')".i "${pp_file}";
    else
      cp "$1" "${pp_file}"
      mv "${modif_file}" "${hash_file}"
    fi
  fi
}

ser_pgi_pp() {
  local pp_file="../${my_iconppdir}/$1"
  local hash_file="../${my_iconppdir}/$1.sha1"
  local modif_file="../${my_iconppdir}/$1.modified"
  if [ -e "${modif_file}" ]; then
    if grep -q '!$ser' $1; then
      ${FC} -I./include ${PP_FFLAGS} -E "$1" > "${pp_file}"
    else
      cp "$1" "${pp_file}"
      mv "${modif_file}" "${hash_file}"
    fi
  fi
}

ser_sb2_pp() {
  local hash_file="$1.sha1"
  local modif_file="$1.modified"
  if [ -e "${modif_file}" ]; then
    ${SERIALBOX2_PPSER} -o "$1" "$1"
    mv "${modif_file}" "${hash_file}"
  fi
}
