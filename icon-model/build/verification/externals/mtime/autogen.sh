#!/bin/sh

# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

set -e

script_dir=`echo "$0" | sed 's@[^/]*$@@'`
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH
cd "$script_dir"

exec ./scripts/reconfigure
