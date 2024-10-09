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

# Checks whether the building generated any git-untracked files.
# Must be run in the root build directory after a successful in-source build.

set -eu
set -o pipefail

icon_dir=$(unset CDPATH; cd "$(dirname "$0")/../../.."; pwd)

exitcode=0

# Check whether the building generated any untracked files in the main repository:
untracked_files=`git -C "${icon_dir}" ls-files --other --exclude-standard` || {
  echo "ERROR: failed to get a list of untracked files in the main repository" >&2
  exit 1
}
if test -n "${untracked_files}"; then
  cat >&2 <<_EOF
ERROR: the building process generated untracked files in the main repository:
(the paths are relative to the source root directory '${icon_dir}')

${untracked_files}

Update '.gitignore' to ignore them.
_EOF
  exitcode=1
fi

# Check whether the building generated any untracked files in the submodules:
untracked_files=`git -C "${icon_dir}" submodule --quiet foreach --recursive 'git ls-files --other --exclude-standard | sed "s|^|$path/|"'` || {
  echo "ERROR: failed to get a list of untracked files in the git submodules" >&2
  exit 1
}
if test -n "${untracked_files}"; then
  cat >&2 <<_EOF
ERROR: the building process generated untracked files in the git submodules:
(the paths are relative to the source root directory '${icon_dir}')

${untracked_files}

Update the respective '.gitignore' files to ignore them.
_EOF
  exitcode=1
fi

exit ${exitcode}
