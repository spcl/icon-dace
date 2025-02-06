# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

pkgconfigdir = $(libdir)/pkgconfig

nodist_pkgconfig_DATA = yac-core.pc

if ENABLE_MCI
nodist_pkgconfig_DATA += yac-mci.pc

if ENABLE_DEPRECATED
nodist_pkgconfig_DATA += yac.pc
endif ENABLE_DEPRECATED

endif ENABLE_MCI

if ENABLE_UTILS
nodist_pkgconfig_DATA += yac-utils.pc
endif ENABLE_UTILS
