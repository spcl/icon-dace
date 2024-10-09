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

# Checks that the input files do not contain any "private" modules from the granules

# The command line arguments can be paths either to regular files or to
# directories. In the latter case, the directories are searched recursively and
# all files with names that match known patterns are considered as input files
# for checking. If no arguments are provided, the script runs the standard
# Buildbot test, i.e. checks files in ICON source directories.

# The known patterns are (newline-separated list of name patterns):
known_patterns='
*.F90
*.f90
*.inc
*.incf
'

# ICON source directories (newline-separated list of paths relative to the root
# source directory):
icon_directories='
src
support
externals
'

# ICON ignored patterns (newline-separated list of path patterns relative to the
# root source directory):
icon_ignored_patterns='
src/granules/microphysics_1mom_schemes/microphysics_1mom_schemes.f90
src/granules/microphysics_1mom_schemes/gscp_data.f90
src/granules/microphysics_1mom_schemes/gscp_kessler.f90
src/granules/microphysics_1mom_schemes/gscp_cloudice.f90
src/granules/microphysics_1mom_schemes/gscp_graupel.f90
src/granules/microphysics_1mom_schemes/gscp_ice.f90
'

# Number of parallel jobs:
job_num=8

set -eu
set -o pipefail

quote_value ()
{
  eval "qv_x=\$${2-${1}}"
  case ${qv_x} in
    *\'*) qv_x=`echo "${qv_x}" | sed "s/'/'\\\\\\\\''/g"` ;;
  esac
  eval "${1}=\"'\${qv_x}'\""
  unset qv_x
}

flatten_and_quote_values ()
{
  eval "faqv_input=\$${2-${1}}; ${1}=''"
  while read faqv_x; do
    test -n "${faqv_x}" || continue
    quote_value faqv_x
    eval "${1}=\"\$${1}${faqv_x} \""
  done <<_EOF
${faqv_input}
_EOF
  unset faqv_input faqv_x
}

check_exist()
{
  for ce_input in "${@}"; do
    if test ! -e "${ce_input}"; then
      echo "ERROR: '${ce_input}' does not exist" >&2
      exit 2
    fi
  done
  unset ce_input
}

list_files()
{
  lf_issue_warn=${1}; shift
  lf_dirs=''
  for lf_input in "${@}"; do
    if test -d "${lf_input}"; then
      quote_value lf_input
      lf_dirs="${lf_dirs} ${lf_input}"
    elif test -f "${lf_input}"; then
      echo "${lf_input}"
    elif test x"${issue_warn}" = xyes; then
      echo "WARNING: input argument '${lf_input}' is neither a directory nor a file" >&2
    fi
  done
  if test -n "${lf_dirs}"; then
    lf_find_args="${lf_dirs} -type f -a \\("
    while read lf_x; do
      test -n "${lf_x}" || continue
      quote_value lf_x
      lf_find_args="${lf_find_args} -name ${lf_x} -o"
    done <<_EOF
${known_patterns}
_EOF
    lf_find_args="${lf_find_args} -false \\) -a \\! \\("
    eval "set dummy ${ignored_patterns}; shift"
    for lf_x in "${@}"; do
      quote_value lf_x
      lf_find_args="${lf_find_args} -path ${lf_x} -o"
    done
    lf_find_args="${lf_find_args} -false \\)"
    eval "find${lf_find_args}"
    unset lf_find_args
  fi
  unset lf_issue_warn lf_dirs lf_input
}

ignored_patterns=''
if test ${#} -eq 0; then
  icon_prefix="$(unset CDPATH; cd "$(dirname "${0}")/../.."; pwd)/"
  quote_value icon_prefix
  flatten_and_quote_values icon_ignored_patterns
  eval "set dummy ${icon_ignored_patterns}; shift"
  for pattern in "${@}"; do
    quote_value pattern
    ignored_patterns="${ignored_patterns} ${icon_prefix}${pattern}"
  done
  flatten_and_quote_values icon_directories
  eval "set dummy ${icon_directories}; shift"
  icon_directories=''
  for dir in "${@}"; do
    quote_value dir
    icon_directories="${icon_directories} ${icon_prefix}${dir}"
  done
  eval "set dummy ${icon_directories}; shift"
fi

check_exist "${@}"

exitcode=0

regexp='gscp_(graupel|cloudice|data|kessler|ice)'
quote_value regexp

list_files no "${@}" | xargs -P ${job_num} -I{} -- ${SHELL-${BASH}} -c "grep --color='auto' -HniP ${regexp} '{}' >&2; test \${?} -eq 1" || {
  {
    echo "ERROR: some of the input files contain private modules from granules"
    echo "       Make sure you only use the public interface of each granule"
  } >&2
  exitcode=1
}

exit ${exitcode}
