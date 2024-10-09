#!/bin/sh

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: CC0-1.0

set -e

script_dir=`echo "$0" | sed 's@[^/]*$@@'`
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH
cd "$script_dir"

exec autoreconf -fvi
