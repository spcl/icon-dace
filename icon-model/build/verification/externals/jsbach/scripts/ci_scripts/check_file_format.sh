#!/usr/bin/env bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#----------------------------------------------------------------------
# Checks that the input files do not contain:
#   1) byte sequences that are not valid UTF-8;
#   2) all UTF-8 characters are also ASCII characters.
#   3) trailing blanks
#   4) too long lines (for fortran code)
# The command line arguments can be paths either to regular files or to
# directories. In the latter case, the directories are searched recursively and
# all files with names that match known patterns are considered as input files
# for checking. If no arguments are provided, the script runs the standard
# Buildbot test, i.e. checks files in ICON source directories.

check_invalid_UTF8=true
check_non_ASCII=true
check_trailing_blanks=true
check_too_long_lines=false

# The known patterns are (space-separated list of single-quoted patterns):
known_patterns="'*f90' '*.F90' '*.inc' '*.incf' '*.h' '*.c' '*.cu'"

# ICON source directories (space-separated list of single-quoted paths relative
# to the root repo directory):
icon_directories="'src' 'scripts'"

# Number of parallel jobs:
job_num=8

set -eu
set -o pipefail

list_files()
{
  issue_warn=$1; shift
  for input in "$@"; do
    if test -d "${input}"; then
      eval "find "${input}" -type f -a \( $(
        eval "set dummy ${known_patterns}; shift"
        for pattern in "$@"; do
          echo -n "-name '${pattern}' -o "
        done
        echo "-false") \)"
    elif test -f "${input}"; then
      echo "${input}"
    elif test x"$issue_warn" = xyes; then
      echo "WARNING: input argument '${input}' is neither a directory nor a file" >&2
    fi
  done
}

if test $# -eq 0; then
  icon_dir=$(unset CDPATH; cd "$(dirname "$0")/../.."; pwd)
  eval "set dummy $(
    eval "set dummy ${icon_directories}; shift"
    for dir in "$@"; do
      echo -n "'${icon_dir}/${dir}' "
    done); shift"
fi

exitcode=0
if [[ ${check_invalid_UTF8} == true ]]; then
  list_files yes "$@" | xargs -P ${job_num} -I{} -- /bin/bash -c 'LC_CTYPE=en_US.UTF-8 grep --color="auto" -Hnaxv ".*" {} >&2; test $? -eq 1' || {
    echo "ERROR: input files contain invalid UTF-8 byte sequences (see above)" >&2
    exitcode=1
  }
fi
if [[ ${check_non_ASCII} == true ]]; then
  list_files no "$@" | xargs -P ${job_num} -I{} -- /bin/bash -c 'grep --color="auto" -HnP "[^\x00-\x7F]" {} >&2; test $? -eq 1' || {
    echo "ERROR: input files contain non-ASCII characters (see above)" >&2
    exitcode=1
  }
fi
if [[ ${check_trailing_blanks} == true ]]; then
  list_files no "$@" | xargs -P ${job_num} -I{} -- /bin/bash -c 'grep --color="auto" -HnP " $" {} >&2; test $? -eq 1' || {
    echo "ERROR: input files contain trailing blanks (see above)" >&2
    exitcode=1
  }
fi
if [[ ${check_too_long_lines} == true ]]; then
  list_files no "$@" | xargs -P ${job_num} -I{} -- /bin/bash -c 'grep --color="auto" -Hn ".\{133,\}" {} >&2; test $? -eq 1' || {
    echo "ERROR: input files contain lines with more then 132 characters (see above)" >&2
    exitcode=1
  }
fi

exit ${exitcode}
