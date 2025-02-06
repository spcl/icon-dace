# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

nodist_include_HEADERS = yac_config.h

SUBDIRS = core utils

if ENABLE_MCI
SUBDIRS += mci

if ENABLE_DEPRECATED
lib_LIBRARIES = libyac.a
# Keep 'make dist' working:
nodist_libyac_a_SOURCES = libyac.c

manifests = core/manifest.txt mci/manifest.txt

libyac.a: $(manifests)
	$(AM_V_at)-rm -f $@
	$(AM_V_AR)$(AR) $(ARFLAGS) $@ `cat $(manifests)`
	$(AM_V_at)$(RANLIB) $@
endif ENABLE_DEPRECATED

endif ENABLE_MCI

SUBDIRS += pkgconfig
