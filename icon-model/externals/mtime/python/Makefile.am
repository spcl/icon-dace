# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
install_tmpdir = .install

pkgpython_PYTHON = mtime/__init__.py
nodist_pkgpython_PYTHON = $(install_tmpdir)/_mtime.py
pkgpyexec_LTLIBRARIES = mtime/__mtime.la

mtime___mtime_la_SOURCES =
mtime___mtime_la_LDFLAGS = -shared -avoid-version -module
mtime___mtime_la_LIBADD = $(top_builddir)/src/libmtime_c.la

$(install_tmpdir)/_mtime.py: mtime/_mtime.py
	$(MKDIR_P) $(@D) && { \
	  $(SED) -n '1,/^# SET-SOFILE-BEGIN-TAG$$/p' $<; \
	  echo '__SOFILE = "$(DESTDIR)$(pkgpyexecdir)/__mtime$(acx_modext)"'; \
	  $(SED) -n '/^# SET-SOFILE-END-TAG$$/,$$p' $<; \
	}>$@

MOSTLYCLEANFILES = $(nodist_pkgpython_PYTHON)

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)
