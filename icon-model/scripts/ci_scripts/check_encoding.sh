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

#/bin/bash

# Checks that the input files do not contain:
#   1) byte sequences that are not valid UTF-8;
#   2) all UTF-8 characters are also ASCII characters.
# The command line arguments can be paths either to regular files or to
# directories. In the latter case, the directories are searched recursively and
# all files with names that match known patterns are considered as input files
# for checking. If no arguments are provided, the script runs the standard
# Buildbot test, i.e. checks files in ICON source directories.

# The known patterns are (space-separated list of single-quoted patterns):
known_patterns="'*.f90' '*.F90' '*.inc' '*.incf' '*.h' '*.c' '*.cu'"

# ICON source directories (space-separated list of single-quoted paths relative
# to the root repo directory):
icon_directories="'src' 'support'"

# Number of parallel jobs:
job_num=8

set -eu
set -o pipefail

check_exist()
{
  for input in "$@"; do
    if test ! -e "${input}"; then
      echo "ERROR: '${input}' does not exist" >&2
      exit 2
    fi
  done
}

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
      echo "WARNING: input argument '${input}' is neither a directory not a file" >&2
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

check_exist "$@"

exitcode=0

list_files yes "$@" | xargs -P ${job_num} -I{} -- ${SHELL-$BASH} -c 'LC_CTYPE=en_US.UTF-8 grep --color="auto" -Hnaxv ".*" {} >&2; test $? -eq 1' || {
  echo "ERROR: input files contain invalid UTF-8 byte sequences (see above)" >&2
  exitcode=1
}

list_files no "$@" | xargs -P ${job_num} -I{} -- ${SHELL-$BASH} -c 'grep --color="auto" -HnP "[^\x00-\x7F]" {} >&2; test $? -eq 1' || {
  echo "ERROR: input files contain non-ASCII characters (see above)" >&2
  exitcode=1
}

exit ${exitcode}
