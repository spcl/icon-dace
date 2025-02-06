# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

try:
    from ._yac import *
except Exception as e:
    print("Your PYTHONPATH probably points to the source directory of yac/python instead of the build directory.")
    raise
